//! Training driver (SPEC.md §10, WP10; `TRAINING.md`).
//!
//! ```text
//! cargo run --release -p ddg-points-as-tori --features train,cuda --example train -- <command> [options]
//!
//! run      the curriculum (training/curriculum.toml) on the authors' data; resumes from --out
//! bench    step timing, memory and run-to-run determinism on one shape (GPU and/or CPU)
//! compare  the same short run on the CPU and the GPU (same data, seed, init): loss curves
//! export   the weights of a checkpoint directory (or a phase file) as a safetensors file
//!
//! options (defaults in brackets):
//!   --data DIR            training data [/var/data/discrete-geometry/PointsAsTori/PAT_training_data/PAT_training_data]
//!   --out DIR             checkpoints, log.jsonl [runs_tmp/pat-train]
//!   --device D            cuda:N | cpu [cuda:0]
//!   --threads T           host threads [all but 4 with a GPU, all on the CPU]
//!   --seed S              [42]
//!   --sign R              reference | geometric [reference]
//!   --eigen R             robust | reference [robust]
//!   --schedule S          code | cosine [code]
//!   --raw-normals         keep the stored norm-2 normals (parity only)
//!   --reference-padding   the reference's padding (parity only)
//!   --max-steps N         stop after N steps in all
//!   --checkpoint-every N  [2000]
//!   --no-validate         skip the validation at epoch ends
//!   --validate-every N    validate at an epoch end once N steps passed since the last one [0: every epoch]
//!   --no-resume           start fresh even if --out holds a checkpoint
//!   --network-batch B     neighborhoods per network launch [256 on the GPU, 4 on the CPU]
//!   --pair-chunk P        padded pairs per loss chunk [1048576 on the GPU, 65536 on the CPU]
//!   --phases a,b          only these phases (names), e.g. pretrain
//!   --file F --shape I --steps N --reps R    (bench, compare)
//!   --from DIR|FILE --to FILE                 (export)
//! ```

use ddg_points_as_tori::fit::{EigenRule, SignRule};
use ddg_points_as_tori::train::data::{DataOptions, load_file};
use ddg_points_as_tori::train::net::NetParams;
use ddg_points_as_tori::train::optim::AdamTrainer;
use ddg_points_as_tori::train::step::{StepOptions, on_gpu, plan_slice, prepare_shape, step_planned};
use ddg_points_as_tori::train::tfit::Ctx;
use ddg_points_as_tori::train::{Schedule, StepLog, TrainOptions};
use ddg_points_as_tori::{Model, PatError};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

const DATA: &str = "/var/data/discrete-geometry/PointsAsTori/PAT_training_data/PAT_training_data";

#[derive(Debug, Clone)]
struct Args {
    cmd: String,
    data: PathBuf,
    out: PathBuf,
    device: String,
    threads: Option<usize>,
    seed: u64,
    sign: SignRule,
    eigen: EigenRule,
    schedule: Schedule,
    raw_normals: bool,
    reference_padding: bool,
    max_steps: Option<u64>,
    checkpoint_every: u64,
    validate: bool,
    validate_every: u64,
    resume: bool,
    network_batch: Option<usize>,
    pair_chunk: Option<usize>,
    phases: Option<Vec<String>>,
    file: String,
    shape: usize,
    steps: usize,
    reps: usize,
    from: Option<PathBuf>,
    to: Option<PathBuf>,
}

fn bad(what: String) -> PatError {
    PatError::InvalidOption { what }
}

fn parse() -> Result<Args, PatError> {
    let mut it = std::env::args().skip(1);
    let cmd = it.next().ok_or_else(|| bad("usage: train <run|bench|compare|export> [options] (see the source header)".into()))?;
    let mut a = Args {
        cmd,
        data: PathBuf::from(DATA),
        out: PathBuf::from("runs_tmp/pat-train"),
        device: "cuda:0".into(),
        threads: None,
        seed: 42,
        sign: SignRule::Reference,
        eigen: EigenRule::Robust,
        schedule: Schedule::Code,
        raw_normals: false,
        reference_padding: false,
        max_steps: None,
        checkpoint_every: 2000,
        validate: true,
        validate_every: 0,
        resume: true,
        network_batch: None,
        pair_chunk: None,
        phases: None,
        file: "nice/blobs_train_mode0_0.npz".into(),
        shape: 0,
        steps: 8,
        reps: 2,
        from: None,
        to: None,
    };
    let args: Vec<String> = it.collect();
    let mut i = 0;
    let val = |i: &mut usize| -> Result<String, PatError> {
        *i += 1;
        args.get(*i).cloned().ok_or_else(|| bad(format!("{} needs a value", args[*i - 1])))
    };
    let num = |s: String, what: &str| -> Result<u64, PatError> { s.parse::<u64>().map_err(|_| bad(format!("{what}: {s:?} is not a number"))) };
    while i < args.len() {
        match args[i].as_str() {
            "--data" => a.data = val(&mut i)?.into(),
            "--out" => a.out = val(&mut i)?.into(),
            "--device" => a.device = val(&mut i)?,
            "--threads" => a.threads = Some(num(val(&mut i)?, "--threads")? as usize),
            "--seed" => a.seed = num(val(&mut i)?, "--seed")?,
            "--sign" => {
                a.sign = match val(&mut i)?.as_str() {
                    "reference" => SignRule::Reference,
                    "geometric" => SignRule::Geometric,
                    s => return Err(bad(format!("--sign {s} (reference | geometric)"))),
                }
            }
            "--eigen" => {
                a.eigen = match val(&mut i)?.as_str() {
                    "robust" => EigenRule::Robust,
                    "reference" => EigenRule::Reference,
                    s => return Err(bad(format!("--eigen {s} (robust | reference)"))),
                }
            }
            "--schedule" => {
                a.schedule = match val(&mut i)?.as_str() {
                    "code" => Schedule::Code,
                    "cosine" => Schedule::Cosine,
                    s => return Err(bad(format!("--schedule {s} (code | cosine)"))),
                }
            }
            "--raw-normals" => a.raw_normals = true,
            "--reference-padding" => a.reference_padding = true,
            "--max-steps" => a.max_steps = Some(num(val(&mut i)?, "--max-steps")?),
            "--checkpoint-every" => a.checkpoint_every = num(val(&mut i)?, "--checkpoint-every")?.max(1),
            "--no-validate" => a.validate = false,
            "--validate-every" => a.validate_every = num(val(&mut i)?, "--validate-every")?,
            "--no-resume" => a.resume = false,
            "--network-batch" => a.network_batch = Some(num(val(&mut i)?, "--network-batch")?.max(1) as usize),
            "--pair-chunk" => a.pair_chunk = Some(num(val(&mut i)?, "--pair-chunk")?.max(1) as usize),
            "--phases" => a.phases = Some(val(&mut i)?.split(',').map(str::to_string).collect()),
            "--file" => a.file = val(&mut i)?,
            "--shape" => a.shape = num(val(&mut i)?, "--shape")? as usize,
            "--steps" => a.steps = num(val(&mut i)?, "--steps")? as usize,
            "--reps" => a.reps = num(val(&mut i)?, "--reps")?.max(1) as usize,
            "--from" => a.from = Some(val(&mut i)?.into()),
            "--to" => a.to = Some(val(&mut i)?.into()),
            s => return Err(bad(format!("unknown option {s}"))),
        }
        i += 1;
    }
    Ok(a)
}

fn ctx_for(device: &str) -> Result<Ctx, PatError> {
    if device == "cpu" {
        return Ok(Ctx { dev: ddg_points_as_tori::network::burn::cpu_device(), dtype: burn::tensor::DType::F32 });
    }
    let ord = device.strip_prefix("cuda:").and_then(|s| s.parse::<usize>().ok()).ok_or_else(|| bad(format!("--device {device} (cuda:N | cpu)")))?;
    cuda(ord)
}

#[cfg(feature = "cuda")]
fn cuda(ord: usize) -> Result<Ctx, PatError> {
    ddg_points_as_tori::train::gpu::cuda_ctx(ord)
}

#[cfg(not(feature = "cuda"))]
fn cuda(_: usize) -> Result<Ctx, PatError> {
    Err(PatError::Gpu { reason: "built without the feature cuda".into() })
}

fn step_options(a: &Args, gpu: bool) -> StepOptions {
    let mut o = StepOptions::default();
    o.loss.fit.sign_rule = a.sign;
    o.loss.fit.eigen_rule = a.eigen;
    o.reference_padding = a.reference_padding;
    o.network_chunk = a.network_batch.unwrap_or(if gpu { 256 } else { 4 });
    o.loss.pair_chunk = a.pair_chunk.unwrap_or(if gpu { 1 << 20 } else { 1 << 16 });
    o
}

fn train_options(a: &Args, gpu: bool) -> Result<TrainOptions, PatError> {
    let mut o = TrainOptions::reference(a.out.clone())?;
    if let Some(keep) = &a.phases {
        o.curriculum.phases.retain(|p| keep.contains(&p.name));
        if o.curriculum.phases.is_empty() {
            return Err(bad(format!("--phases {keep:?}: no such phase")));
        }
    }
    o.step = step_options(a, gpu);
    o.data = DataOptions { raw_normals: a.raw_normals, ..Default::default() };
    o.schedule = a.schedule;
    o.seed = a.seed;
    o.checkpoint_every = a.checkpoint_every;
    o.max_steps = a.max_steps;
    o.resume = a.resume;
    o.validate = a.validate;
    o.validate_every_steps = a.validate_every;
    o.threads = a.threads;
    Ok(o)
}

fn host_threads(a: &Args, gpu: bool) -> Option<usize> {
    a.threads.or_else(|| {
        let n = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
        // leave cores for burn's dispatch threads (NOTES-GPU.md: three threads at 80-90 %)
        gpu.then(|| n.saturating_sub(4).max(1))
    })
}

fn fmt_log(r: &StepLog, t0: std::time::Instant, total: Option<u64>, start: u64) -> String {
    let el = t0.elapsed().as_secs_f64();
    let done = r.step.saturating_sub(start).max(1);
    let eta = total.map(|t| (t.saturating_sub(r.step)) as f64 * el / done as f64 / 3600.0);
    match r.validation {
        Some((d, e)) => format!("[{:>9.0} s] {} epoch {} validation: distance {d:.5}, eikonal {e:.5}", el, r.phase_name, r.epoch),
        None => format!(
            "[{:>9.0} s] step {:>7} {} e{} {} loss {:.5} (d {:.5}, e {:.5}) lr {:.1e} pairs {} rows {}/{} t {:.3}+{:.3}+{:.3}+{:.3} s{}{}",
            el,
            r.step,
            r.phase_name,
            r.epoch,
            r.shape,
            r.loss,
            r.distance,
            r.eikonal,
            r.lr,
            r.pairs,
            r.rows.0,
            r.rows.1,
            r.seconds[2],
            r.seconds[3],
            r.seconds[4],
            r.seconds[5],
            if r.skipped { " SKIPPED" } else { "" },
            eta.map(|h| format!(" eta {h:.1} h")).unwrap_or_default()
        ),
    }
}

fn run(a: &Args) -> Result<(), PatError> {
    let ctx = ctx_for(&a.device)?;
    let gpu = on_gpu(&ctx);
    let mut o = train_options(a, gpu)?;
    o.threads = host_threads(a, gpu);
    let pool = ddg_points_as_tori::pool::make_pool(o.threads)?;
    let total: u64 = 497_858;
    eprintln!("train: device {} ({}), host threads {:?}, out {}, options {:?}", a.device, if gpu { "GPU" } else { "CPU" }, o.threads, a.out.display(), o.step);
    let t0 = std::time::Instant::now();
    let mut first: Option<u64> = None;
    let mut last_print = std::time::Instant::now() - std::time::Duration::from_secs(3600);
    let (model, report) = ddg_points_as_tori::train::run::train_on(&a.data, None, &o, &ctx, &pool, &mut |r: &StepLog| {
        let s = *first.get_or_insert(r.step.saturating_sub(1));
        if r.validation.is_some() || r.skipped || last_print.elapsed().as_secs_f64() >= 30.0 || r.step <= s + 5 {
            eprintln!("{}", fmt_log(r, t0, (a.phases.is_none()).then_some(total), s));
            last_print = std::time::Instant::now();
        }
    })?;
    eprintln!("train: done, {report:?}");
    let w = NetParams::from_model(&model, &ctx)
        .to_safetensors(&HashMap::from([("source".to_string(), "ddg-points-as-tori training (WP10), final".to_string())]))?;
    std::fs::write(a.out.join("final.safetensors"), w).map_err(|e| bad(format!("final.safetensors: {e}")))?;
    Ok(())
}

fn hash_f32(v: &[f32]) -> u64 {
    v.iter().fold(0xcbf29ce484222325u64, |h, x| (h ^ x.to_bits() as u64).wrapping_mul(0x100000001b3))
}

fn grads_hash(g: &[ddg_points_as_tori::train::net::PT]) -> u64 {
    hash_f32(&g.iter().flat_map(|t| t.to_f32()).collect::<Vec<f32>>())
}

fn params_hash(p: &NetParams) -> u64 {
    hash_f32(&p.items.iter().flat_map(|t| t.value.to_f32()).collect::<Vec<f32>>())
}

fn nvsmi() -> String {
    std::process::Command::new("nvidia-smi")
        .args(["--query-gpu=memory.used,utilization.gpu", "--format=csv,noheader"])
        .output()
        .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        .unwrap_or_default()
}

/// Step timing, memory and run-to-run determinism on one shape.
fn bench(a: &Args) -> Result<(), PatError> {
    let ctx = ctx_for(&a.device)?;
    let gpu = on_gpu(&ctx);
    let o = step_options(a, gpu);
    let pool = ddg_points_as_tori::pool::make_pool(host_threads(a, gpu))?;
    let (shapes, _) = pool.install(|| load_file(&a.data.join(&a.file), &DataOptions::default()))?;
    let shape = shapes.get(a.shape).ok_or_else(|| bad(format!("{} has {} shapes", a.file, shapes.len())))?;
    let t = std::time::Instant::now();
    let prep = pool.install(|| prepare_shape(shape, &o, (a.seed, 0, a.shape as u64)))?;
    let t_prep = t.elapsed().as_secs_f64();
    let batch = if a.file.contains("pretraining") { 32 } else { 1024 };
    let t = std::time::Instant::now();
    let plans: Vec<_> = (0..shape.rows.min(shape.n.max(1)))
        .step_by(batch)
        .map(|lo| pool.install(|| plan_slice(&prep, lo, lo + batch, &o, gpu)))
        .collect::<Result<Vec<_>, _>>()?;
    let plans: Vec<_> = plans.into_iter().flatten().collect();
    let t_plan = t.elapsed().as_secs_f64() / plans.len().max(1) as f64;
    {
        // the parts of the plan: blend constants (kd-tree) and the GPU pair layout
        let t = std::time::Instant::now();
        let sp = pool.install(|| plan_slice(&prep, 0, batch, &o, false))?;
        let t_bp = t.elapsed().as_secs_f64();
        #[cfg(feature = "cuda")]
        if let Some(sp) = &sp {
            let t = std::time::Instant::now();
            let l = ddg_points_as_tori::train::gpu::pair_layout(&sp.plan, prep.valid.len(), &o.loss);
            eprintln!(
                "plan parts: blend constants {t_bp:.3} s, GPU layout {:.3} s ({} query chunks, {} torus chunks, {} padded pairs)",
                t.elapsed().as_secs_f64(),
                l.chunks.len(),
                l.tchunks.len(),
                l.m_pad
            );
        }
        let _ = (sp, t_bp);
    }
    eprintln!(
        "bench {} shape {} ({} points, {} valid): queries + truth {t_prep:.3} s, plan {t_plan:.3} s per slice ({} slices, {} pairs in the first), device {}, batch {}, pair chunk {}",
        a.file,
        a.shape,
        shape.rows,
        shape.n,
        plans.len(),
        plans.first().map(|p| p.plan.pairs()).unwrap_or(0),
        a.device,
        o.network_chunk,
        o.loss.pair_chunk
    );
    // the authors' weights (a trained regime) and the initialization
    for (name, init) in [("authors", NetParams::from_model(&Model::reference()?, &ctx)), ("init", NetParams::init(&Default::default(), a.seed, &ctx))] {
        let mut hashes = Vec::new();
        for rep in 0..a.reps {
            let mut params = init.clone();
            let mut adam = AdamTrainer::new(&params, Default::default());
            let mut line = Vec::new();
            let t_all = std::time::Instant::now();
            for s in 0..a.steps {
                let sp = &plans[s % plans.len()];
                let t = std::time::Instant::now();
                let out =
                    if gpu { step_planned(&params, &ctx, &prep, sp, &o, true)? } else { pool.install(|| step_planned(&params, &ctx, &prep, sp, &o, true))? };
                let g = out.grads.ok_or_else(|| bad("no gradient".into()))?;
                let gh = grads_hash(&g);
                let tu = std::time::Instant::now();
                params = adam.update(&params, &g, 1e-5)?;
                let ph = params_hash(&params);
                let t_upd = tu.elapsed().as_secs_f64();
                eprintln!(
                    "  {name} rep {rep} step {s}: loss {:.7} (d {:.7} e {:.7}) bits {:016x} grads {gh:016x} params {ph:016x}; fwd {:.3} loss {:.3} bwd {:.3} upd {t_upd:.3} = {:.3} s; rows {}/{}; GPU {}",
                    out.value.loss,
                    out.value.distance,
                    out.value.eikonal,
                    out.value.loss.to_bits(),
                    out.seconds[1],
                    out.seconds[2],
                    out.seconds[3],
                    t.elapsed().as_secs_f64(),
                    out.rows_forward,
                    out.rows_backward,
                    if gpu { nvsmi() } else { String::new() }
                );
                line.push((out.value.loss.to_bits(), gh, ph));
            }
            eprintln!("  {name} rep {rep}: {} steps in {:.2} s", a.steps, t_all.elapsed().as_secs_f64());
            hashes.push(line);
        }
        let same = hashes.windows(2).all(|w| w[0] == w[1]);
        eprintln!("{name}: {} repetitions of {} steps bitwise identical: {same}", a.reps, a.steps);
    }
    Ok(())
}

fn read_log(p: &Path) -> Vec<(u64, f64, f64, f64)> {
    let s = std::fs::read_to_string(p).unwrap_or_default();
    let num = |l: &str, k: &str| -> Option<f64> {
        let pat = format!("\"{k}\":");
        let at = l.find(&pat)? + pat.len();
        let r = &l[at..];
        r[..r.find([',', '}']).unwrap_or(r.len())].parse().ok()
    };
    s.lines()
        .filter(|l| !l.contains("val_distance"))
        .filter_map(|l| Some((num(l, "step")? as u64, num(l, "loss")?, num(l, "distance")?, num(l, "eikonal")?)))
        .collect()
}

/// The same short run (phases given by --phases, default pretraining) on the CPU and the GPU.
/// One step from the same parameters on the CPU and the GPU: loss and gradient differences.
fn compare_step(a: &Args, file: &str, batch: usize) -> Result<(), PatError> {
    let (cpu, gpu) = (ctx_for("cpu")?, ctx_for(&a.device)?);
    let pool = ddg_points_as_tori::pool::make_pool(host_threads(a, true))?;
    let (shapes, _) = pool.install(|| load_file(&a.data.join(file), &DataOptions::default()))?;
    let shape = shapes.get(a.shape).ok_or_else(|| bad(format!("{file} has {} shapes", shapes.len())))?;
    let (oc, og) = (step_options(a, false), step_options(a, true));
    let prep = pool.install(|| prepare_shape(shape, &oc, (a.seed, 0, a.shape as u64)))?;
    let sp = pool.install(|| plan_slice(&prep, 0, batch, &og, true))?.ok_or_else(|| bad("empty slice".into()))?;
    for (name, model) in [("authors", Model::reference()?), ("init", NetParams::init(&Default::default(), a.seed, &cpu).to_model())] {
        let t = std::time::Instant::now();
        let c = pool.install(|| step_planned(&NetParams::from_model(&model, &cpu), &cpu, &prep, &sp, &oc, true))?;
        let tc = t.elapsed().as_secs_f64();
        let t = std::time::Instant::now();
        let g = step_planned(&NetParams::from_model(&model, &gpu), &gpu, &prep, &sp, &og, true)?;
        let tg = t.elapsed().as_secs_f64();
        let flat =
            |o: &ddg_points_as_tori::train::StepOut| -> Vec<f64> { o.grads.as_ref().map(|g| g.iter().flat_map(|t| t.to_f64()).collect()).unwrap_or_default() };
        let (gc, gg) = (flat(&c), flat(&g));
        let norm = |v: &[f64]| v.iter().map(|x| x * x).sum::<f64>().sqrt();
        let diff: Vec<f64> = gc.iter().zip(&gg).map(|(x, y)| x - y).collect();
        let gmax = gc.iter().fold(0f64, |m, x| m.max(x.abs()));
        let dmax = diff.iter().fold(0f64, |m, x| m.max(x.abs()));
        let cos = gc.iter().zip(&gg).map(|(x, y)| x * y).sum::<f64>() / (norm(&gc) * norm(&gg)).max(1e-300);
        eprintln!(
            "step {file} shape {} rows 0..{batch}, {name}: loss CPU {:.9} GPU {:.9} (rel {:.2e}; distance rel {:.2e}, eikonal rel {:.2e}); gradient |Δ|/|g| {:.2e}, max|Δ|/max|g| {:.2e}, cosine {:.8}; CPU {tc:.2} s, GPU {tg:.2} s",
            a.shape,
            c.value.loss,
            g.value.loss,
            (c.value.loss - g.value.loss).abs() / c.value.loss.abs(),
            (c.value.distance - g.value.distance).abs() / c.value.distance.abs(),
            (c.value.eikonal - g.value.eikonal).abs() / c.value.eikonal.abs(),
            norm(&diff) / norm(&gc).max(1e-300),
            dmax / gmax.max(1e-300),
            cos
        );
    }
    Ok(())
}

fn compare(a: &Args) -> Result<(), PatError> {
    compare_step(a, "pretraining_blob.npz", 32)?;
    compare_step(a, &a.file, 1024)?;
    let mut a = a.clone();
    if a.phases.is_none() {
        a.phases = Some(vec!["pretrain".into()]);
    }
    a.max_steps = Some(a.steps as u64);
    a.validate = false;
    a.resume = false;
    let mut logs = Vec::new();
    for dev in ["cuda:0", "cpu"] {
        let ctx = ctx_for(dev)?;
        let gpu = on_gpu(&ctx);
        let mut b = a.clone();
        b.out = a.out.join(if gpu { "gpu" } else { "cpu" });
        let _ = std::fs::remove_file(b.out.join("log.jsonl"));
        let _ = std::fs::remove_file(b.out.join("checkpoint.safetensors"));
        let mut o = train_options(&b, gpu)?;
        o.threads = host_threads(&b, gpu);
        let pool = ddg_points_as_tori::pool::make_pool(o.threads)?;
        let t0 = std::time::Instant::now();
        let (_, rep) =
            ddg_points_as_tori::train::run::train_on(&b.data, None, &o, &ctx, &pool, &mut |r: &StepLog| eprintln!("{dev}: {}", fmt_log(r, t0, None, 0)))?;
        eprintln!("{dev}: {} steps in {:.1} s ({:?})", rep.steps, rep.seconds, rep);
        logs.push(read_log(&b.out.join("log.jsonl")));
    }
    let (g, c) = (&logs[0], &logs[1]);
    let mut worst = 0f64;
    eprintln!("step  loss GPU      loss CPU      rel. diff");
    for (x, y) in g.iter().zip(c) {
        let r = (x.1 - y.1).abs() / y.1.abs().max(1e-30);
        worst = worst.max(r);
        eprintln!("{:>4}  {:.9}  {:.9}  {r:.2e}  (d {:.3e} e {:.3e})", x.0, x.1, y.1, (x.2 - y.2).abs(), (x.3 - y.3).abs());
    }
    eprintln!("compare: {} steps, max relative loss difference GPU vs CPU {worst:.3e}", g.len().min(c.len()));
    Ok(())
}

/// Export the weights of a checkpoint directory or a phase file.
fn export(a: &Args) -> Result<(), PatError> {
    let from = a.from.clone().unwrap_or_else(|| a.out.clone());
    let to = a.to.clone().ok_or_else(|| bad("export needs --to FILE".into()))?;
    let model = if from.is_dir() {
        ddg_points_as_tori::train::run::checkpoint_model(&from)?
    } else {
        let b = std::fs::read(&from).map_err(|e| bad(format!("{}: {e}", from.display())))?;
        Model::from_safetensors(&b)?
    };
    let ctx = ctx_for("cpu")?;
    let mut meta = HashMap::new();
    meta.insert("source".to_string(), format!("ddg-points-as-tori retraining (WP10), exported from {}", from.display()));
    if let Some(st) = model.metadata().get("train_state") {
        meta.insert("train_state".to_string(), st.clone());
    }
    meta.insert(
        "license".to_string(),
        "Trained by ddg-points-as-tori with a port of points-as-tori's training code (MIT License, Copyright (c) 2026 Nicole Feng) on the authors' training data.".into(),
    );
    let w = NetParams::from_model(&model, &ctx).to_safetensors(&meta)?;
    if let Some(d) = to.parent() {
        std::fs::create_dir_all(d).map_err(|e| bad(format!("{}: {e}", d.display())))?;
    }
    std::fs::write(&to, &w).map_err(|e| bad(format!("{}: {e}", to.display())))?;
    // round trip
    let back = Model::from_safetensors(&w)?;
    eprintln!("export: {} -> {} ({} bytes, {} parameters)", from.display(), to.display(), w.len(), NetParams::from_model(&back, &ctx).count());
    Ok(())
}

fn main() {
    let r = parse().and_then(|a| match a.cmd.as_str() {
        "run" => run(&a),
        "bench" => bench(&a),
        "compare" => compare(&a),
        "export" => export(&a),
        c => Err(bad(format!("unknown command {c} (run | bench | compare | export)"))),
    });
    if let Err(e) = r {
        eprintln!("train: {e}");
        std::process::exit(1);
    }
}
