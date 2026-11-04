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
use super::step::{StepOptions, evaluate_shape, prepare_shape, step};
use super::tfit::Ctx;
use crate::{Model, ModelConfig, PatError};
use burn::tensor::DType;
use rand_chacha::rand_core::RngCore;
use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

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
            if batch_rows == 0 || !(lr > 0.0) {
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
fn epoch_order(counts: &[(usize, usize)], seed: u64, phase: usize, epoch: u64, per_buffer: usize) -> (Vec<(usize, usize)>, Vec<(usize, usize)>) {
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
pub fn validate(params: &NetParams, ctx: &Ctx, file: &Path, data: &DataOptions, o: &StepOptions, batch: usize, seed: u64, epoch: u64) -> Result<(f64, f64), PatError> {
    let (shapes, _) = load_file(file, data)?;
    let (mut sd, mut se, mut c) = (0.0, 0.0, 0usize);
    for (i, s) in shapes.iter().enumerate() {
        let Ok(prep) = prepare_shape(s, o, (seed, epoch, 10_000 + i as u64)) else { continue };
        let v = evaluate_shape(params, ctx, &prep, batch, o)?;
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

/// [`train`] on a given burn device and dtype (WP10 passes a CUDA device); the parallel host work
/// runs on `pool`.
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
    let mut gstep = 0u64;
    let start_phase = resume.as_ref().map(|r| r.2.phase).unwrap_or(0);
    // weights of the phase before the resumed one
    if start_phase > 0 && resume.is_none() {
        return Err(PatError::InvalidOption { what: "cannot start at a later phase without a checkpoint".into() });
    }
    let mut stop = false;
    for (pi, ph) in opts.curriculum.phases.iter().enumerate().skip(start_phase) {
        let files = phase_files(data_dir, ph)?;
        let counts = shape_counts(&files)?;
        let steps_per_epoch: u64 = counts.iter().map(|&(s, p)| (s * p.div_ceil(ph.batch_rows)) as u64).sum();
        let phase_steps = steps_per_epoch * ph.epochs;
        let (mut adam, mut cur) = match resume.take() {
            Some((p, a, c)) => {
                params = p;
                gstep = c.step;
                (a, c)
            }
            None => {
                if ph.from_scratch && pi > 0 {
                    params = match init {
                        Some(m) => NetParams::from_model(m, ctx),
                        None => NetParams::init(&cfg, opts.seed, ctx),
                    };
                }
                (AdamTrainer::new(&params, opts.optimizer), Cursor { phase: pi, step: gstep, ..Default::default() })
            }
        };
        let mut since_ckpt = 0u64;
        for epoch in cur.epoch..ph.epochs {
            let (seq, bufs) = epoch_order(&counts, opts.seed, pi, epoch, opts.files_in_buffer);
            let skip = if epoch == cur.epoch { cur.done } else { 0 };
            for &(blo, bhi) in &bufs {
                if bhi <= skip {
                    continue;
                }
                // load the buffer's files
                let mut want: Vec<usize> = seq[blo..bhi].iter().map(|x| x.0).collect();
                want.sort_unstable();
                want.dedup();
                let mut loaded = HashMap::new();
                for f in want {
                    let (shapes, dr) = pool.install(|| load_file(&files[f], &opts.data))?;
                    report.data.merge(&dr);
                    loaded.insert(f, shapes);
                }
                for (pos, &(f, si)) in seq.iter().enumerate().take(bhi).skip(blo.max(skip)) {
                    let shape = &loaded[&f][si];
                    let sname = format!("{}:{}", shape.source.0, si);
                    let tp = std::time::Instant::now();
                    let gid = 100_000 * f as u64 + si as u64;
                    let prep = match pool.install(|| prepare_shape(shape, &opts.step, (opts.seed, 1_000_000 * pi as u64 + epoch, gid))) {
                        Ok(p) => p,
                        Err(_) => {
                            report.skipped_shapes.push(sname);
                            continue;
                        }
                    };
                    report.rows_without_faces += prep.queries.rows_without_faces;
                    let mut t_prep = tp.elapsed().as_secs_f64();
                    for lo in (0..shape.rows).step_by(ph.batch_rows) {
                        let hi = lo + ph.batch_rows;
                        let mut rec = StepLog { phase: pi, phase_name: ph.name.clone(), epoch, shape: sname.clone(), row: lo, ..Default::default() };
                        let lr = learning_rate(opts.schedule, adam.count, ph.lr, ph.warmup, phase_steps);
                        rec.lr = lr;
                        if lo >= shape.n {
                            if !opts.keep_empty_slices {
                                report.skipped_empty += 1;
                                continue;
                            }
                            let zero: Vec<_> = params.items.iter().map(|p| p.value.scale(0.0)).collect();
                            params = adam.update(&params, &zero, lr)?;
                        } else {
                            let out = pool.install(|| step(&params, ctx, &prep, lo, hi, &opts.step, true))?;
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
                        t_prep = 0.0;
                        gstep += 1;
                        since_ckpt += 1;
                        report.steps += (!rec.skipped) as u64;
                        rec.step = gstep;
                        let _ = writeln!(log, "{}", rec.json());
                        progress(&rec);
                        if opts.max_steps.is_some_and(|m| gstep >= m) {
                            stop = true;
                            break;
                        }
                    }
                    cur = Cursor { phase: pi, epoch, done: pos + 1, step: gstep, ..cur };
                    if since_ckpt >= opts.checkpoint_every || stop {
                        save_checkpoint(&opts.checkpoint_dir, &params, &adam, &cur, "")?;
                        since_ckpt = 0;
                    }
                    if stop {
                        break;
                    }
                }
                if stop {
                    break;
                }
            }
            if stop {
                break;
            }
            // end of epoch
            cur = Cursor { phase: pi, epoch: epoch + 1, done: 0, step: gstep, ..cur };
            save_checkpoint(&opts.checkpoint_dir, &params, &adam, &cur, "")?;
            since_ckpt = 0;
            if opts.validate
                && let Some(v) = &opts.curriculum.validation
            {
                let val = pool.install(|| validate(&params, ctx, &data_dir.join(v), &opts.data, &opts.step, 1024, opts.seed, epoch))?;
                report.validation = Some(val);
                let rec = StepLog { phase: pi, phase_name: ph.name.clone(), epoch, step: gstep, validation: Some(val), ..Default::default() };
                let _ = writeln!(log, "{}", rec.json());
                progress(&rec);
            }
        }
        if stop {
            break;
        }
        // end of phase: keep its weights, the next phase starts at its first epoch
        let w = params.to_safetensors(&HashMap::from([("source".to_string(), format!("ddg-points-as-tori training, phase {}", ph.name))]))?;
        write_atomic(&opts.checkpoint_dir.join(format!("phase-{}.safetensors", ph.name)), &w)?;
        let next = Cursor { phase: pi + 1, epoch: 0, done: 0, step: gstep, ..Default::default() };
        save_checkpoint(&opts.checkpoint_dir, &params, &AdamTrainer::new(&params, opts.optimizer), &next, "")?;
    }
    report.seconds = t0.elapsed().as_secs_f64();
    Ok((params.to_model(), report))
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
