"""Run the extra-normaliser and bootstrap jobs for the six real matrices under a memory
budget (estimated peak GB per job; jobs start when the running estimate + the job's estimate
stays below BUDGET_GB). Logs go to logs/. Usage: python3 run_all.py [--skip-done]"""
import os, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
BUDGET_GB = 22
ENV = dict(os.environ, OMP_NUM_THREADS="4", OPENBLAS_NUM_THREADS="4", MKL_NUM_THREADS="4")
# (script, matrix index, tag, estimated peak GB)  -- sizes from nnz: HCCp 79M, HCCe 142M, COADp 136M,
# COADe 238M, LUADp 196M, LUADe 488M
JOBS = [("benchmark_bootstrap.py", 5, "LUAD_enhancer", 17),
        ("benchmark_extra_normalisers.py", 5, "LUAD_enhancer", 12),
        ("benchmark_bootstrap.py", 3, "COAD_enhancer", 9),
        ("benchmark_extra_normalisers.py", 3, "COAD_enhancer", 6),
        ("benchmark_bootstrap.py", 4, "LUAD_proximal", 7),
        ("benchmark_extra_normalisers.py", 4, "LUAD_proximal", 5),
        ("benchmark_bootstrap.py", 2, "COAD_proximal", 4),
        ("benchmark_extra_normalisers.py", 2, "COAD_proximal", 4),
        ("benchmark_bootstrap.py", 1, "HCC_enhancer", 4),
        ("benchmark_extra_normalisers.py", 1, "HCC_enhancer", 4),
        ("benchmark_bootstrap.py", 0, "HCC_proximal", 3),
        ("benchmark_extra_normalisers.py", 0, "HCC_proximal", 3)]


def is_running(script, idx):
    r = subprocess.run(["pgrep", "-f", f"{script} --matrix {idx}$"], capture_output=True, text=True)
    return bool(r.stdout.strip())


def done(script, tag):
    p = HERE / "parts" / (("boot_" if "bootstrap" in script else "extra_") + tag + (".npz" if "bootstrap" in script else ".csv"))
    return p.exists()


def main():
    skip = "--skip-done" in sys.argv
    pending = list(JOBS)
    running, running_ext = [], []
    while pending or running or running_ext:
        running_ext = [j for j in running_ext if is_running(j[0], j[1])]
        for proc, job in list(running):
            if proc.poll() is not None:
                running.remove((proc, job))
                print(f"[{time.strftime('%H:%M:%S')}] finished {job[0]} {job[2]} rc={proc.returncode}", flush=True)
        used = sum(j[3] for _, j in running) + sum(j[3] for j in running_ext)
        for job in list(pending):
            if skip and done(job[0], job[2]):
                pending.remove(job); print(f"skip {job[0]} {job[2]} (done)", flush=True); continue
            if skip and is_running(job[0], job[1]):
                pending.remove(job); used += job[3]; running_ext.append(job)
                print(f"already running {job[0]} {job[2]} (budget used {used} GB)", flush=True); continue
            if used + job[3] <= BUDGET_GB:
                logf = HERE / "logs" / f"{'boot' if 'bootstrap' in job[0] else 'extra'}_{job[2]}.log"
                proc = subprocess.Popen([sys.executable, str(HERE / job[0]), "--matrix", str(job[1])],
                                        stdout=open(logf, "w"), stderr=subprocess.STDOUT, env=ENV, cwd=HERE)
                running.append((proc, job)); pending.remove(job); used += job[3]
                print(f"[{time.strftime('%H:%M:%S')}] started {job[0]} {job[2]} (budget used {used} GB)", flush=True)
        time.sleep(15)
    print("all jobs finished", flush=True)


if __name__ == "__main__":
    main()
