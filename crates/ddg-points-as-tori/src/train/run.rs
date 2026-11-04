//! The training loop (SPEC.md §10): the curriculum (`training/curriculum.toml`), data order (D22),
//! steps, checkpoints (safetensors), logs (JSON lines) and resumption.
//!
//! **Order.** Files of a phase are listed per pattern and sorted by name (the reference uses
//! unsorted `glob`, N§9.2). Each epoch permutes the files, takes them in buffers of
//! [`TrainOptions::files_in_buffer`] and permutes the shapes of each buffer (ChaCha8 keyed by seed,
//! phase and epoch). With every file in one buffer (pretraining, phase 0) this is the reference's
//! global shuffle; for the large phases it bounds memory (phase 4 would otherwise hold ≈ 16 GB).
//!
//! **Steps.** As the reference: per shape, `⌈P/batch_rows⌉` steps over consecutive row slices, each
//! running the network on the cloud. Slices without valid rows are skipped (counted) unless
//! [`TrainOptions::keep_empty_slices`]; a step with a non-finite loss or gradient is skipped and
//! counted (SPEC.md §10.2); under `reference_padding` non-finite gradient entries are zeroed
//! instead, as the reference.
//!
//! **Checkpoints** are written at shape boundaries (every [`TrainOptions::checkpoint_every`] steps,
//! at the end of every epoch and phase): `checkpoint.safetensors` (weights, readable by
//! [`crate::Model::from_safetensors`], with the cursor in its metadata) and
//! `checkpoint-adam.safetensors` (moments). A resumed run continues at the next shape and is
//! identical to an uninterrupted one (queries are keyed by (seed, epoch, shape)).

use super::data::{DataOptions, DataReport, load_file};
use super::net::NetParams;
use super::npz::Npz;
use super::optim::{AdamOptions, AdamTrainer, Schedule, learning_rate};
use super::queries::row_rng;
use super::data::TrainShape;
use super::step::{PreparedShape, SlicePlan, StepOptions, StepOut, on_gpu, plan_slice, prepare_shape, step_planned};
use super::tfit::Ctx;
use crate::{Model, ModelConfig, PatError};
use burn::tensor::DType;
use rand_chacha::rand_core::RngCore;
use std::collections::HashMap;
use std::io::Write;
use rayon::prelude::*;
use std::path::{Path, PathBuf};
use std::sync::{Arc, mpsc};

/// One phase of the curriculum.
#[derive(Debug, Clone, PartialEq)]
pub struct Phase {
    /// Name (`pretrain`, `phase0`, …).
    pub name: String,
    /// Groups of file patterns relative to the data directory (`*` in the file name only); each
    /// pattern's matches are sorted by name.
    pub files: Vec<String>,
    /// Interleave the patterns' lists, truncated to the shortest (the reference's `zip` in phase 1).
    pub interleave: bool,
    /// Epochs.
    pub epochs: u64,
    /// Peak learning rate.
    pub lr: f64,
    /// Linear warmup steps.
    pub warmup: u64,
    /// Query rows per step (1024; 32 in pretraining).
    pub batch_rows: usize,
    /// Start from fresh parameters (`true` for the pretraining phase) or from the previous phase.
    pub from_scratch: bool,
}

/// The curriculum.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Curriculum {
    /// Phases in order.
    pub phases: Vec<Phase>,
    /// Validation file relative to the data directory (optional).
    pub validation: Option<String>,
}

fn terr(what: String) -> PatError {
    PatError::Format { reason: format!("curriculum: {what}") }
}

impl Curriculum {
    /// Parse `training/curriculum.toml`.
    pub fn from_toml(s: &str) -> Result<Curriculum, PatError> {
        let t: toml::Table = s.parse().map_err(|e| terr(format!("{e}")))?;
        let validation = t.get("validation").and_then(|v| v.as_str()).map(str::to_string);
        let phases = t.get("phase").and_then(|v| v.as_array()).ok_or_else(|| terr("no [[phase]] tables".into()))?;
        let mut out = Vec::new();
        for (i, p) in phases.iter().enumerate() {
            let p = p.as_table().ok_or_else(|| terr(format!("phase {i} is not a table")))?;
            let get = |k: &str| p.get(k).ok_or_else(|| terr(format!("phase {i} lacks {k}")));
            let int = |k: &str| -> Result<u64, PatError> {
                get(k)?.as_integer().filter(|&v| v >= 0).map(|v| v as u64).ok_or_else(|| terr(format!("phase {i}: {k} is not a non-negative integer")))
            };
            let lr = get("lr")?.as_float().ok_or_else(|| terr(format!("phase {i}: lr is not a float")))?;
            let files = get("files")?
                .as_array()
                .ok_or_else(|| terr(format!("phase {i}: files is not an array")))?
                .iter()
                .map(|f| f.as_str().map(str::to_string).ok_or_else(|| terr(format!("phase {i}: a file pattern is not a string"))))
                .collect::<Result<Vec<_>, _>>()?;
            let init = p.get("init").and_then(|v| v.as_str()).unwrap_or("previous");
            if init != "scratch" && init != "previous" {
                return Err(terr(format!("phase {i}: init = {init:?} (scratch or previous)")));
            }
            let batch_rows = int("batch_rows")? as usize;
            if batch_rows == 0 || lr.is_nan() || lr <= 0.0 {
                return Err(terr(format!("phase {i}: batch_rows and lr must be positive")));
            }
            out.push(Phase {
                name: get("name")?.as_str().ok_or_else(|| terr(format!("phase {i}: name is not a string")))?.to_string(),
                files,
                interleave: p.get("interleave").and_then(|v| v.as_bool()).unwrap_or(false),
                epochs: int("epochs")?,
                lr,
                warmup: int("warmup")?,
                batch_rows,
                from_scratch: init == "scratch",
            });
        }
        Ok(Curriculum { phases: out, validation })
    }

    /// The curriculum committed with the crate (`training/curriculum.toml`).
    pub fn reference() -> Result<Curriculum, PatError> {
        Curriculum::from_toml(include_str!("../../training/curriculum.toml"))
    }
}

/// `*` matches any run of characters (no other metacharacters).
fn wildcard(pat: &str, s: &str) -> bool {
    let parts: Vec<&str> = pat.split('*').collect();
    if parts.len() == 1 {
        return pat == s;
    }
    let (first, last) = (parts[0], parts[parts.len() - 1]);
    if !s.starts_with(first) || !s[first.len()..].ends_with(last) || s.len() < first.len() + last.len() {
        return false;
    }
    let mut rest = &s[first.len()..s.len() - last.len()];
    for p in &parts[1..parts.len() - 1] {
        match rest.find(p) {
            Some(i) => rest = &rest[i + p.len()..],
            None => return false,
        }
    }
    true
}

/// Expand `dir/name` patterns (`*` wildcards in the file name only), sorted.
pub fn expand(data_dir: &Path, pattern: &str) -> Result<Vec<PathBuf>, PatError> {
    let p = data_dir.join(pattern);
    let (dir, name) = (p.parent().unwrap_or(data_dir).to_path_buf(), p.file_name().map(|f| f.to_string_lossy().to_string()).unwrap_or_default());
    if !name.contains('*') {
        return if p.is_file() { Ok(vec![p]) } else { Err(PatError::Format { reason: format!("training file {} not found", p.display()) }) };
    }
    let mut out: Vec<PathBuf> = std::fs::read_dir(&dir)
        .map_err(|e| PatError::Format { reason: format!("{}: {e}", dir.display()) })?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|f| f.file_name().map(|n| n.to_string_lossy().to_string()).is_some_and(|n| wildcard(&name, &n)))
        .collect();
    out.sort();
    if out.is_empty() {
        return Err(PatError::Format { reason: format!("no training file matches {}", p.display()) });
    }
    Ok(out)
}

/// The files of a phase.
pub fn phase_files(data_dir: &Path, ph: &Phase) -> Result<Vec<PathBuf>, PatError> {
    let lists = ph.files.iter().map(|f| expand(data_dir, f)).collect::<Result<Vec<_>, _>>()?;
    if ph.interleave {
        let m = lists.iter().map(|l| l.len()).min().unwrap_or(0);
        Ok((0..m).flat_map(|i| lists.iter().map(move |l| l[i].clone())).collect())
    } else {
        Ok(lists.concat())
    }
}

/// Options of [`train`].
#[derive(Debug, Clone, PartialEq)]
pub struct TrainOptions {
    /// The phases.
    pub curriculum: Curriculum,
    /// Steps (loss, queries, padding, chunks).
    pub step: StepOptions,
    /// Data loading (raw normals).
    pub data: DataOptions,
    /// Adam/AdamW.
    pub optimizer: AdamOptions,
    /// Schedule after the warmup (D21).
    pub schedule: Schedule,
    /// Seed (data order, queries, initialization).
    pub seed: u64,
    /// Checkpoints and `log.jsonl`.
    pub checkpoint_dir: PathBuf,
    /// Steps between checkpoints (at the next shape boundary).
    pub checkpoint_every: u64,
    /// Stop after this many steps in all (short runs).
    pub max_steps: Option<u64>,
    /// Step Adam on slices without valid rows, as the reference (zero gradient).
    pub keep_empty_slices: bool,
    /// Files per shuffle buffer.
    pub files_in_buffer: usize,
    /// Continue from `checkpoint_dir` if it holds a checkpoint.
    pub resume: bool,
    /// Evaluate the curriculum's validation file at the end of every epoch.
    pub validate: bool,
    /// Validate at an epoch end only when at least this many steps were taken since the last
    /// validation (always at the last epoch of a phase); 0: every epoch.
    pub validate_every_steps: u64,
    /// Threads (`None`: all cores).
    pub threads: Option<usize>,
}

impl TrainOptions {
    /// The reference's setup with the committed curriculum.
    pub fn reference(checkpoint_dir: PathBuf) -> Result<TrainOptions, PatError> {
        Ok(TrainOptions {
            curriculum: Curriculum::reference()?,
            step: StepOptions::default(),
            data: DataOptions::default(),
            optimizer: AdamOptions::default(),
            schedule: Schedule::Code,
            seed: 42,
            checkpoint_dir,
            checkpoint_every: 1000,
            max_steps: None,
            keep_empty_slices: false,
            files_in_buffer: 4,
            resume: true,
            validate: true,
            validate_every_steps: 0,
            threads: None,
        })
    }
}

/// One step's record (also a line of `log.jsonl`).
#[derive(Debug, Clone, PartialEq, Default)]
pub struct StepLog {
    /// Phase index and name.
    pub phase: usize,
    /// Phase name.
    pub phase_name: String,
    /// Epoch within the phase.
    pub epoch: u64,
    /// Global step (all phases).
    pub step: u64,
    /// Shape (`file:index`).
    pub shape: String,
    /// First query row of the slice.
    pub row: usize,
    /// Learning rate used.
    pub lr: f64,
    /// Loss, distance and eikonal terms before the update.
    pub loss: f64,
    /// Distance term.
    pub distance: f64,
    /// Eikonal term.
    pub eikonal: f64,
    /// Queries.
    pub queries: usize,
    /// Pairs in the blend.
    pub pairs: usize,
    /// Network rows forward / backward.
    pub rows: (usize, usize),
    /// Step seconds: queries+GT (shared per shape, on its first step), plan, forward, loss,
    /// backward, update.
    pub seconds: [f64; 6],
    /// Skipped (non-finite).
    pub skipped: bool,
    /// Validation (distance, eikonal), at epoch ends.
    pub validation: Option<(f64, f64)>,
}

impl StepLog {
    fn json(&self) -> String {
        let v = self.validation.map(|(d, e)| format!(",\"val_distance\":{d:e},\"val_eikonal\":{e:e}")).unwrap_or_default();
        format!(
            "{{\"phase\":{},\"phase_name\":\"{}\",\"epoch\":{},\"step\":{},\"shape\":\"{}\",\"row\":{},\"lr\":{:e},\"loss\":{:e},\"distance\":{:e},\"eikonal\":{:e},\"queries\":{},\"pairs\":{},\"rows_fwd\":{},\"rows_bwd\":{},\"t_prep\":{:.4},\"t_plan\":{:.4},\"t_fwd\":{:.4},\"t_loss\":{:.4},\"t_bwd\":{:.4},\"t_upd\":{:.4},\"skipped\":{}{v}}}",
            self.phase, self.phase_name, self.epoch, self.step, self.shape, self.row, self.lr, self.loss, self.distance, self.eikonal,
            self.queries, self.pairs, self.rows.0, self.rows.1, self.seconds[0], self.seconds[1], self.seconds[2], self.seconds[3], self.seconds[4], self.seconds[5], self.skipped
        )
    }
}

/// What training did.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct TrainReport {
    /// Steps taken (updates).
    pub steps: u64,
    /// Steps skipped for a non-finite loss or gradient.
    pub skipped_nonfinite: u64,
    /// Slices without valid rows that were skipped.
    pub skipped_empty: u64,
    /// Shapes that could not be used (no face of positive area), as `file:index`.
    pub skipped_shapes: Vec<String>,
    /// Rows whose regularized box met no face.
    pub rows_without_faces: usize,
    /// Data loading.
    pub data: DataReport,
    /// Last validation (distance, eikonal).
    pub validation: Option<(f64, f64)>,
    /// Resumed from a checkpoint at this global step.
    pub resumed_at: Option<u64>,
    /// Wall time.
    pub seconds: f64,
}

/// The cursor of a checkpoint.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
struct Cursor {
    phase: usize,
    epoch: u64,
    done: usize,
    step: u64,
    adam_count: u64,
    adam_time: u64,
}

fn shape_counts(files: &[PathBuf]) -> Result<Vec<(usize, usize)>, PatError> {
    files
        .iter()
        .map(|f| {
            let mut z = Npz::open(f)?;
            let (s, _) = z.read_bool("masks")?;
            if s.len() != 2 {
                return Err(PatError::Format { reason: format!("{}: masks has shape {s:?}", f.display()) });
            }
            Ok((s[0], s[1]))
        })
        .collect()
}

fn permutation(n: usize, key: (u64, u64, u64)) -> Vec<usize> {
    let mut r = row_rng(key.0, key.1, key.2, u64::MAX);
    let mut v: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        let j = (r.next_u64() % (i as u64 + 1)) as usize;
        v.swap(i, j);
    }
    v
}

/// The epoch's sequence of `(file, shape)` and its buffers (ranges into the sequence).
/// `(file, shape)` pairs, or `[lo, hi)` ranges.
type Pairs = Vec<(usize, usize)>;

fn epoch_order(counts: &[(usize, usize)], seed: u64, phase: usize, epoch: u64, per_buffer: usize) -> (Pairs, Pairs) {
    let fo = permutation(counts.len(), (seed, 1_000_000 * phase as u64 + epoch, 0xf11e));
    let mut seq = Vec::new();
    let mut bufs = Vec::new();
    for (b, chunk) in fo.chunks(per_buffer.max(1)).enumerate() {
        let items: Vec<(usize, usize)> = chunk.iter().flat_map(|&f| (0..counts[f].0).map(move |s| (f, s))).collect();
        let perm = permutation(items.len(), (seed, 1_000_000 * phase as u64 + epoch, 0xb0f0 + b as u64));
        let lo = seq.len();
        seq.extend(perm.into_iter().map(|i| items[i]));
        bufs.push((lo, seq.len()));
    }
    (seq, bufs)
}

fn write_atomic(path: &Path, bytes: &[u8]) -> Result<(), PatError> {
    let tmp = path.with_extension("tmp");
    std::fs::write(&tmp, bytes).map_err(|e| PatError::Format { reason: format!("{}: {e}", tmp.display()) })?;
    std::fs::rename(&tmp, path).map_err(|e| PatError::Format { reason: format!("{}: {e}", path.display()) })
}

fn save_checkpoint(dir: &Path, params: &NetParams, adam: &AdamTrainer, cur: &Cursor, extra: &str) -> Result<(), PatError> {
    let (named, time) = adam.to_named(params);
    let mut meta = HashMap::new();
    meta.insert(
        "train_state".to_string(),
        format!(
            "{{\"phase\": {}, \"epoch\": {}, \"done\": {}, \"step\": {}, \"adam_count\": {}, \"adam_time\": {}}}",
            cur.phase, cur.epoch, cur.done, cur.step, adam.count, time
        ),
    );
    meta.insert("source".into(), format!("ddg-points-as-tori training (WP7){extra}"));
    meta.insert(
        "license".into(),
        "Trained by ddg-points-as-tori with a port of points-as-tori's training code (MIT License, Copyright (c) 2026 Nicole Feng).".into(),
    );
    let w = params.to_safetensors(&meta)?;
    let bytes: Vec<(String, Vec<usize>, Vec<u8>)> = named.into_iter().map(|(n, s, v)| (n, s, v.iter().flat_map(|x| x.to_le_bytes()).collect())).collect();
    let views: Vec<(String, safetensors::tensor::TensorView<'_>)> = bytes
        .iter()
        .map(|(n, s, b)| safetensors::tensor::TensorView::new(safetensors::Dtype::F32, s.clone(), b).map(|v| (n.clone(), v)).map_err(|e| PatError::Weights { reason: format!("{n}: {e}") }))
        .collect::<Result<_, _>>()?;
    let a = safetensors::serialize(views, None::<HashMap<String, String>>).map_err(|e| PatError::Weights { reason: format!("{e}") })?;
    write_atomic(&dir.join("checkpoint-adam.safetensors"), &a)?;
    write_atomic(&dir.join("checkpoint.safetensors"), &w)
}

fn json_u64(s: &str, k: &str) -> Option<u64> {
    let pat = format!("\"{k}\":");
    let at = s.find(&pat)? + pat.len();
    let r = s[at..].trim_start();
    r[..r.find(|c: char| !c.is_ascii_digit()).unwrap_or(r.len())].parse().ok()
}

fn load_checkpoint(dir: &Path, ctx: &Ctx, opts: &AdamOptions) -> Result<Option<(NetParams, AdamTrainer, Cursor)>, PatError> {
    let wp = dir.join("checkpoint.safetensors");
    if !wp.is_file() {
        return Ok(None);
    }
    let wb = std::fs::read(&wp).map_err(|e| PatError::Format { reason: format!("{}: {e}", wp.display()) })?;
    let model = Model::from_safetensors(&wb)?;
    let st = model.metadata().get("train_state").cloned().ok_or_else(|| PatError::Format { reason: format!("{}: no train_state", wp.display()) })?;
    let g = |k: &str| json_u64(&st, k).ok_or_else(|| PatError::Format { reason: format!("{}: train_state lacks {k}", wp.display()) });
    let cur = Cursor { phase: g("phase")? as usize, epoch: g("epoch")?, done: g("done")? as usize, step: g("step")?, adam_count: g("adam_count")?, adam_time: g("adam_time")? };
    let params = NetParams::from_model(&model, ctx);
    let ab = std::fs::read(dir.join("checkpoint-adam.safetensors")).map_err(|e| PatError::Format { reason: format!("adam state: {e}") })?;
    let st = safetensors::SafeTensors::deserialize(&ab).map_err(|e| PatError::Format { reason: format!("adam state: {e}") })?;
    let named: HashMap<String, Vec<f32>> = st
        .tensors()
        .into_iter()
        .map(|(n, v)| (n, v.data().chunks_exact(4).map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]])).collect()))
        .collect();
    let adam = AdamTrainer::from_named(&params, *opts, &named, cur.adam_time, cur.adam_count, ctx)?;
    Ok(Some((params, adam, cur)))
}

/// Mean validation (distance, eikonal) over the shapes of `file`, as the reference: per shape the
/// mean over its batches, then the mean over shapes; queries keyed by `(seed, epoch, 10000 + i)`.
/// Runs on the current rayon pool.
#[allow(clippy::too_many_arguments)]
pub fn validate(params: &NetParams, ctx: &Ctx, file: &Path, data: &DataOptions, o: &StepOptions, batch: usize, seed: u64, epoch: u64) -> Result<(f64, f64), PatError> {
    let (shapes, _) = load_file(file, data)?;
    validate_shapes(params, ctx, &shapes, o, batch, seed, epoch, &|f| f())
}

/// [`validate`] with the host work on `pool`; on a CUDA device the network runs on the calling
/// thread (one dispatching thread), on the CPU inside `pool`.
#[allow(clippy::too_many_arguments)]
pub fn validate_on(params: &NetParams, ctx: &Ctx, pool: &rayon::ThreadPool, file: &Path, data: &DataOptions, o: &StepOptions, batch: usize, seed: u64, epoch: u64) -> Result<(f64, f64), PatError> {
    let (shapes, _) = pool.install(|| load_file(file, data))?;
    let host = |f: &(dyn Fn() -> Result<(), PatError> + Sync)| pool.install(f);
    if on_gpu(ctx) {
        validate_shapes(params, ctx, &shapes, o, batch, seed, epoch, &host)
    } else {
        pool.install(|| validate_shapes(params, ctx, &shapes, o, batch, seed, epoch, &|f| f()))
    }
}

type HostRunner<'a> = dyn Fn(&(dyn Fn() -> Result<(), PatError> + Sync)) -> Result<(), PatError> + 'a;

#[allow(clippy::too_many_arguments)]
fn validate_shapes(params: &NetParams, ctx: &Ctx, shapes: &[TrainShape], o: &StepOptions, batch: usize, seed: u64, epoch: u64, host: &HostRunner<'_>) -> Result<(f64, f64), PatError> {
    let (mut sd, mut se, mut c) = (0.0, 0.0, 0usize);
    let gpu = on_gpu(ctx);
    for (i, s) in shapes.iter().enumerate() {
        let cell = std::sync::Mutex::new(None);
        host(&|| {
            let p = prepare_shape(s, o, (seed, epoch, 10_000 + i as u64)).ok();
            let mut v = Vec::new();
            if let Some(p) = &p {
                let mut lo = 0;
                while lo < p.rows {
                    if lo < p.n {
                        v.push(plan_slice(p, lo, lo + batch, o, gpu)?);
                    }
                    lo += batch;
                }
            }
            *cell.lock().map_err(|_| PatError::Backend { reason: "poisoned lock".into() })? = Some((p, v));
            Ok(())
        })?;
        let Some((Some(prep), plans)) = cell.into_inner().map_err(|_| PatError::Backend { reason: "poisoned lock".into() })? else { continue };
        let mut v = Vec::new();
        for sp in plans.iter().flatten() {
            v.push(step_planned(params, ctx, &prep, sp, o, false)?.value);
        }
        if v.is_empty() {
            continue;
        }
        sd += v.iter().map(|x| x.distance).sum::<f64>() / v.len() as f64;
        se += v.iter().map(|x| x.eikonal).sum::<f64>() / v.len() as f64;
        c += 1;
    }
    let c = c.max(1) as f64;
    Ok((sd / c, se / c))
}

/// Train on the CPU (SPEC.md §10.4 API). `init` replaces the parameters of a `from_scratch` phase
/// (fine-tuning); `progress` sees every step.
pub fn train(data_dir: &Path, init: Option<&Model>, opts: &TrainOptions, progress: &mut dyn FnMut(&StepLog)) -> Result<(Model, TrainReport), PatError> {
    let ctx = Ctx { dev: crate::network::burn::cpu_device(), dtype: DType::F32 };
    let pool = crate::pool::make_pool(opts.threads)?;
    train_on(data_dir, init, opts, &ctx, &pool, progress)
}

/// Messages of the loader thread (whole shuffle buffers) to the preparer.
enum LoadMsg {
    Phase { pi: usize, phase_steps: u64 },
    Buffer { pi: usize, epoch: u64, items: Vec<(usize, usize, usize)>, shapes: HashMap<usize, Vec<TrainShape>>, report: DataReport },
    EpochEnd { pi: usize, epoch: u64 },
    PhaseEnd { pi: usize },
}

/// Events of the preparer thread to the training loop, in data order.
enum Ev {
    Phase { pi: usize, phase_steps: u64 },
    Data(DataReport),
    Skipped(String),
    Shape { rows_without_faces: usize },
    Slice { pi: usize, epoch: u64, pos: usize, name: String, lo: usize, prep: Arc<PreparedShape>, plan: Option<Box<SlicePlan>>, empty: bool, t_prep: f64 },
    ShapeDone { pi: usize, epoch: u64, pos: usize },
    EpochEnd { pi: usize, epoch: u64 },
    PhaseEnd { pi: usize },
}

/// Slices prepared ahead of the training loop (queries, ground truth and plans).
const PREFETCH_SLICES: usize = 4;

/// The loader: per phase and epoch the data order, and every shuffle buffer's files (in parallel on
/// `pool`), one buffer ahead of the preparer.
fn loader(data_dir: &Path, opts: &TrainOptions, start: Cursor, pool: &rayon::ThreadPool, tx: &mpsc::SyncSender<Result<LoadMsg, PatError>>) -> Result<(), PatError> {
    let send = |m: LoadMsg| tx.send(Ok(m)).map_err(|_| PatError::Backend { reason: "stopped".into() });
    for (pi, ph) in opts.curriculum.phases.iter().enumerate().skip(start.phase) {
        let files = phase_files(data_dir, ph)?;
        let counts = shape_counts(&files)?;
        let steps_per_epoch: u64 = counts.iter().map(|&(s, p)| (s * p.div_ceil(ph.batch_rows)) as u64).sum();
        send(LoadMsg::Phase { pi, phase_steps: steps_per_epoch * ph.epochs })?;
        let first_epoch = if pi == start.phase { start.epoch } else { 0 };
        for epoch in first_epoch..ph.epochs {
            let (seq, bufs) = epoch_order(&counts, opts.seed, pi, epoch, opts.files_in_buffer);
            let skip = if pi == start.phase && epoch == start.epoch { start.done } else { 0 };
            for &(blo, bhi) in &bufs {
                if bhi <= skip {
                    continue;
                }
                let mut want: Vec<usize> = seq[blo..bhi].iter().map(|x| x.0).collect();
                want.sort_unstable();
                want.dedup();
                let loaded: Vec<Result<(Vec<TrainShape>, DataReport), PatError>> = pool.install(|| want.par_iter().map(|&f| load_file(&files[f], &opts.data)).collect());
                let mut shapes = HashMap::new();
                let mut report = DataReport::default();
                for (&f, l) in want.iter().zip(loaded) {
                    let (s, dr) = l?;
                    report.merge(&dr);
                    shapes.insert(f, s);
                }
                let items = (blo.max(skip)..bhi).map(|pos| (pos, seq[pos].0, seq[pos].1)).collect();
                send(LoadMsg::Buffer { pi, epoch, items, shapes, report })?;
            }
            send(LoadMsg::EpochEnd { pi, epoch })?;
        }
        send(LoadMsg::PhaseEnd { pi })?;
    }
    Ok(())
}

/// The preparer: queries, ground truth and the host plan of every slice, a few slices ahead.
fn preparer(opts: &TrainOptions, gpu: bool, pool: &rayon::ThreadPool, rx: mpsc::Receiver<Result<LoadMsg, PatError>>, tx: &mpsc::SyncSender<Result<Ev, PatError>>) -> Result<(), PatError> {
    let send = |e: Ev| tx.send(Ok(e)).map_err(|_| PatError::Backend { reason: "stopped".into() });
    let mut batch_rows = 1024;
    for m in rx {
        match m? {
            LoadMsg::Phase { pi, phase_steps } => {
                batch_rows = opts.curriculum.phases[pi].batch_rows;
                send(Ev::Phase { pi, phase_steps })?;
            }
            LoadMsg::EpochEnd { pi, epoch } => send(Ev::EpochEnd { pi, epoch })?,
            LoadMsg::PhaseEnd { pi } => send(Ev::PhaseEnd { pi })?,
            LoadMsg::Buffer { pi, epoch, items, shapes, report } => {
                send(Ev::Data(report))?;
                for (pos, f, si) in items {
                    let shape = &shapes[&f][si];
                    let name = format!("{}:{}", shape.source.0, si);
                    let tp = std::time::Instant::now();
                    let gid = 100_000 * f as u64 + si as u64;
                    let prep = match pool.install(|| prepare_shape(shape, &opts.step, (opts.seed, 1_000_000 * pi as u64 + epoch, gid))) {
                        Ok(p) => Arc::new(p),
                        Err(_) => {
                            send(Ev::Skipped(name))?;
                            send(Ev::ShapeDone { pi, epoch, pos })?;
                            continue;
                        }
                    };
                    send(Ev::Shape { rows_without_faces: prep.queries.rows_without_faces })?;
                    let mut t_prep = tp.elapsed().as_secs_f64();
                    for lo in (0..shape.rows).step_by(batch_rows) {
                        let empty = lo >= shape.n;
                        let plan = if empty { None } else { pool.install(|| plan_slice(&prep, lo, lo + batch_rows, &opts.step, gpu))?.map(Box::new) };
                        send(Ev::Slice { pi, epoch, pos, name: name.clone(), lo, prep: prep.clone(), plan, empty, t_prep })?;
                        t_prep = 0.0;
                    }
                    send(Ev::ShapeDone { pi, epoch, pos })?;
                }
            }
        }
    }
    Ok(())
}

/// [`train`] on a given burn device and dtype (a CUDA context from
/// [`super::gpu::cuda_ctx`] with the feature `cuda`); the parallel host work runs on `pool`.
///
/// The data are prefetched: a loader thread reads the next shuffle buffer's files while the current
/// one trains, and a preparer thread computes the queries, ground truth and blend plans
/// ([`plan_slice`]) up to [`PREFETCH_SLICES`] steps ahead, so the device does not wait for the host
/// work, which does not depend on the parameters. The steps run on the calling thread on a CUDA
/// device (one dispatching thread) and inside `pool` on the CPU. Results do not depend on the
/// prefetching.
pub fn train_on(
    data_dir: &Path,
    init: Option<&Model>,
    opts: &TrainOptions,
    ctx: &Ctx,
    pool: &rayon::ThreadPool,
    progress: &mut dyn FnMut(&StepLog),
) -> Result<(Model, TrainReport), PatError> {
    let t0 = std::time::Instant::now();
    std::fs::create_dir_all(&opts.checkpoint_dir).map_err(|e| PatError::Format { reason: format!("{}: {e}", opts.checkpoint_dir.display()) })?;
    let mut log = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(opts.checkpoint_dir.join("log.jsonl"))
        .map_err(|e| PatError::Format { reason: format!("log.jsonl: {e}") })?;
    let mut report = TrainReport::default();
    let cfg = init.map(|m| *m.config()).unwrap_or_default();
    let mut resume = if opts.resume { load_checkpoint(&opts.checkpoint_dir, ctx, &opts.optimizer)? } else { None };
    if let Some((_, _, c)) = &resume {
        report.resumed_at = Some(c.step);
    }
    let mut params = match init {
        Some(m) => NetParams::from_model(m, ctx),
        None => NetParams::init(&cfg, opts.seed, ctx),
    };
    let start = resume.as_ref().map(|r| r.2).unwrap_or_default();
    if start.phase >= opts.curriculum.phases.len() {
        report.seconds = t0.elapsed().as_secs_f64();
        if let Some((p, _, _)) = resume {
            params = p;
        }
        return Ok((params.to_model(), report));
    }
    let gpu = on_gpu(ctx);
    let mut gstep = start.step;
    let (ltx, lrx) = mpsc::sync_channel::<Result<LoadMsg, PatError>>(1);
    let (etx, erx) = mpsc::sync_channel::<Result<Ev, PatError>>(PREFETCH_SLICES);
    std::thread::scope(|sc| -> Result<(Model, TrainReport), PatError> {
        let etx2 = etx.clone();
        sc.spawn(move || {
            if let Err(e) = loader(data_dir, opts, start, pool, &ltx) {
                let _ = ltx.send(Err(e));
            }
        });
        sc.spawn(move || {
            if let Err(e) = preparer(opts, gpu, pool, lrx, &etx2) {
                let _ = etx2.send(Err(e));
            }
        });
        drop(etx);
        let mut adam = AdamTrainer::new(&params, opts.optimizer);
        let mut cur = start;
        let mut phase_steps = 0u64;
        let mut since_ckpt = 0u64;
        let mut last_val = gstep;
        // a burn-cuda failure (out of memory, a lost device) surfaces as a panic at the next read
        // back: it is caught here and returned as an error (the run resumes from its checkpoint)
        let caught = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| -> Result<(), PatError> {
            for ev in &erx {
                let ev = match ev {
                    Ok(ev) => ev,
                    Err(PatError::Backend { reason }) if reason == "stopped" => break,
                    Err(e) => return Err(e),
                };
                match ev {
                    Ev::Phase { pi, phase_steps: ps } => {
                        phase_steps = ps;
                        let ph = &opts.curriculum.phases[pi];
                        match resume.take() {
                            Some((p, a, c)) => {
                                params = p;
                                gstep = c.step;
                                adam = a;
                                cur = c;
                            }
                            None => {
                                if ph.from_scratch && pi > 0 {
                                    params = match init {
                                        Some(m) => NetParams::from_model(m, ctx),
                                        None => NetParams::init(&cfg, opts.seed, ctx),
                                    };
                                }
                                adam = AdamTrainer::new(&params, opts.optimizer);
                                cur = Cursor { phase: pi, step: gstep, ..Default::default() };
                            }
                        }
                        since_ckpt = 0;
                    }
                    Ev::Data(dr) => report.data.merge(&dr),
                    Ev::Skipped(name) => report.skipped_shapes.push(name),
                    Ev::Shape { rows_without_faces } => report.rows_without_faces += rows_without_faces,
                    Ev::Slice { pi, epoch, pos, name, lo, prep, plan, empty, t_prep } => {
                        let ph = &opts.curriculum.phases[pi];
                        let mut rec = StepLog { phase: pi, phase_name: ph.name.clone(), epoch, shape: name, row: lo, ..Default::default() };
                        let lr = learning_rate(opts.schedule, adam.count, ph.lr, ph.warmup, phase_steps);
                        rec.lr = lr;
                        if empty {
                            if !opts.keep_empty_slices {
                                report.skipped_empty += 1;
                                continue;
                            }
                            let zero: Vec<_> = params.items.iter().map(|p| p.value.scale(0.0)).collect();
                            params = adam.update(&params, &zero, lr)?;
                        } else {
                            let out = match &plan {
                                Some(sp) if gpu => step_planned(&params, ctx, &prep, sp, &opts.step, true)?,
                                Some(sp) => pool.install(|| step_planned(&params, ctx, &prep, sp, &opts.step, true))?,
                                None => StepOut::default(),
                            };
                            rec.loss = out.value.loss;
                            rec.distance = out.value.distance;
                            rec.eikonal = out.value.eikonal;
                            rec.queries = out.value.queries;
                            rec.pairs = out.pairs;
                            rec.rows = (out.rows_forward, out.rows_backward);
                            rec.seconds[1..5].copy_from_slice(&out.seconds);
                            let bad = out.value.nonfinite > 0 || !out.value.loss.is_finite() || (out.nonfinite_grads > 0 && !opts.step.reference_padding);
                            if bad || out.grads.is_none() {
                                report.skipped_nonfinite += 1;
                                rec.skipped = true;
                            } else if let Some(g) = out.grads {
                                let tu = std::time::Instant::now();
                                params = adam.update(&params, &g, lr)?;
                                rec.seconds[5] = tu.elapsed().as_secs_f64();
                            }
                        }
                        rec.seconds[0] = t_prep;
                        gstep += 1;
                        since_ckpt += 1;
                        report.steps += (!rec.skipped) as u64;
                        rec.step = gstep;
                        let _ = writeln!(log, "{}", rec.json());
                        progress(&rec);
                        if opts.max_steps.is_some_and(|m| gstep >= m) {
                            // as an interrupted shape: the checkpoint resumes at the next one
                            cur = Cursor { phase: pi, epoch, done: pos + 1, step: gstep, ..cur };
                            save_checkpoint(&opts.checkpoint_dir, &params, &adam, &cur, "")?;
                            break;
                        }
                    }
                    Ev::ShapeDone { pi, epoch, pos } => {
                        cur = Cursor { phase: pi, epoch, done: pos + 1, step: gstep, ..cur };
                        if since_ckpt >= opts.checkpoint_every {
                            save_checkpoint(&opts.checkpoint_dir, &params, &adam, &cur, "")?;
                            since_ckpt = 0;
                        }
                    }
                    Ev::EpochEnd { pi, epoch } => {
                        let ph = &opts.curriculum.phases[pi];
                        cur = Cursor { phase: pi, epoch: epoch + 1, done: 0, step: gstep, ..cur };
                        save_checkpoint(&opts.checkpoint_dir, &params, &adam, &cur, "")?;
                        since_ckpt = 0;
                        let due = epoch + 1 == ph.epochs || gstep.saturating_sub(last_val) >= opts.validate_every_steps;
                        if opts.validate
                            && due
                            && let Some(v) = &opts.curriculum.validation
                        {
                            last_val = gstep;
                            let val = validate_on(&params, ctx, pool, &data_dir.join(v), &opts.data, &opts.step, 1024, opts.seed, epoch)?;
                            report.validation = Some(val);
                            let rec = StepLog { phase: pi, phase_name: ph.name.clone(), epoch, step: gstep, validation: Some(val), ..Default::default() };
                            let _ = writeln!(log, "{}", rec.json());
                            progress(&rec);
                        }
                    }
                    Ev::PhaseEnd { pi } => {
                        let ph = &opts.curriculum.phases[pi];
                        // end of phase: keep its weights, the next phase starts at its first epoch
                        let w = params.to_safetensors(&HashMap::from([("source".to_string(), format!("ddg-points-as-tori training, phase {}", ph.name))]))?;
                        write_atomic(&opts.checkpoint_dir.join(format!("phase-{}.safetensors", ph.name)), &w)?;
                        let next = Cursor { phase: pi + 1, epoch: 0, done: 0, step: gstep, ..Default::default() };
                        save_checkpoint(&opts.checkpoint_dir, &params, &AdamTrainer::new(&params, opts.optimizer), &next, "")?;
                    }
                }
            }
            Ok(())
        }));
        // stop the producers (they see the closed channel at their next send)
        drop(erx);
        match caught {
            Ok(r) => r?,
            Err(e) => {
                let msg = e.downcast_ref::<String>().cloned().or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string())).unwrap_or_default();
                return Err(PatError::Backend { reason: format!("the training step panicked (device error?): {msg}") });
            }
        }
        report.seconds = t0.elapsed().as_secs_f64();
        Ok((params.to_model(), report))
    })
}

/// The model of a checkpoint directory (its `checkpoint.safetensors`).
pub fn checkpoint_model(dir: &Path) -> Result<Model, PatError> {
    let p = dir.join("checkpoint.safetensors");
    let b = std::fs::read(&p).map_err(|e| PatError::Format { reason: format!("{}: {e}", p.display()) })?;
    Model::from_safetensors(&b)
}

/// Default architecture (for docs).
pub fn default_config() -> ModelConfig {
    ModelConfig::default()
}
