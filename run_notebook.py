"""
Run main_playground_clean.ipynb (or any notebook) top-to-bottom without user input.

- Executes every code cell in order in a fresh Jupyter kernel (same as "Restart & Run All").
- Plots are rendered off-screen (MPLBACKEND=Agg), so no windows pop up and nothing waits for you.
  Figures are still saved to figures/ by the notebook code itself.
- After every cell, the executed notebook (with outputs) is saved, so you can open it and
  see exactly how far the run got, even if the job is killed.
- Writes a report (Markdown + JSON) with the status, runtime and error of every cell.

Usage (run from the MayoMate folder, with the environment activated):

    python run_notebook.py                          # run everything, stop at the first error
    python run_notebook.py --continue-on-error      # run everything, record errors, keep going
    python run_notebook.py --skip 3 43              # skip cells by number (see --list / the report)
    python run_notebook.py --only 2 19-29           # run only these cells (ranges allowed)
    python run_notebook.py --list                   # print cell indices and first lines, then exit
    python run_notebook.py --rerun-all              # reprocess every sample (ignore the sample cache)
    python run_notebook.py --rerun JMT49_CK JMT50_CK  # reprocess only these samples

Only new samples are processed: the slow per-sample steps (import, filters, ploidy, recombination
maps, JT clusters) are cached in outputs/sample_cache/ the first time a sample is run, and loaded from
there on later runs. The first run processes everything; later runs only process samples that were
added (or whose input file changed). Analyses that combine samples (filtering of common mutations,
statistics, figures) always use all samples. Use --rerun-all after changing how samples are processed
in the mayo code. Details: mayo/sample_cache.py.

Cell numbers count from 1 at the top of the notebook, including text (markdown) cells,
i.e. the way you count cells in VS Code. Use --list to see them; the report uses the same numbers.
"""

import argparse
import json
import os
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError, CellTimeoutError, DeadKernelError

ANSI_ESCAPE = __import__("re").compile(r"\x1b\[[0-9;]*m")


def parse_indices(values):
    """Turn ["2", "5-7"] into {2, 5, 6, 7}."""
    out = set()
    for v in values or []:
        for part in str(v).split(","):
            part = part.strip()
            if not part:
                continue
            if "-" in part:
                a, b = part.split("-", 1)
                out.update(range(int(a), int(b) + 1))
            else:
                out.add(int(part))
    return out


def first_line(cell):
    for line in cell.source.splitlines():
        if line.strip():
            return line.strip()[:100]
    return ""


def fmt_seconds(s):
    s = int(round(s))
    h, rem = divmod(s, 3600)
    m, sec = divmod(rem, 60)
    return f"{h}h{m:02d}m{sec:02d}s" if h else (f"{m}m{sec:02d}s" if m else f"{sec}s")


def error_text(cell):
    """Return 'ErrorName: message' from a cell's outputs, if it has an error output."""
    for out in cell.get("outputs", []):
        if out.get("output_type") == "error":
            return f"{out.get('ename')}: {ANSI_ESCAPE.sub('', str(out.get('evalue')))}"
    return None


def archive_logs(log_dir: Path, dest: Path, since: float) -> list:
    """Copy log files in log_dir that were written after `since` (a time.time() value) into dest.
    The logs in outputs/logs are overwritten on every run; these copies are kept per run."""
    import shutil
    copied = []
    if not Path(log_dir).is_dir():
        return copied
    for f in sorted(Path(log_dir).iterdir()):
        if f.is_file() and f.stat().st_mtime >= since - 1:
            dest.mkdir(parents=True, exist_ok=True)
            shutil.copy2(f, dest / f.name)
            copied.append(f.name)
    return copied


class MemoryWatcher:
    """Samples the kernel's memory (incl. child processes) every 2 s; tracks the peak per cell and overall."""

    def __init__(self, pid):
        import threading
        try:
            import psutil
            self.proc = psutil.Process(pid)
        except Exception:
            self.proc = None
        self.cell_peak = 0.0
        self.overall_peak = 0.0
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)
        if self.proc is not None:
            self._thread.start()

    def _rss_gb(self):
        try:
            procs = [self.proc] + self.proc.children(recursive=True)
            return sum(p.memory_info().rss for p in procs) / 1024 ** 3
        except Exception:
            return 0.0

    def _run(self):
        while not self._stop.wait(2):
            gb = self._rss_gb()
            self.cell_peak = max(self.cell_peak, gb)
            self.overall_peak = max(self.overall_peak, gb)

    def start_cell(self):
        self.cell_peak = self._rss_gb()

    def stop(self):
        self._stop.set()


def write_report(results, args, started, finished, report_md, report_json, executed_path, peak_gb=None):
    n_ok = sum(r["status"] == "ok" for r in results)
    n_err = sum(r["status"] == "error" for r in results)
    n_skip = sum(r["status"] == "skipped" for r in results)
    n_pending = sum(r["status"] == "not run" for r in results)
    total = sum(r["seconds"] for r in results)

    with open(report_json, "w", encoding="utf-8") as fh:
        json.dump(
            {
                "notebook": str(args.notebook),
                "executed_notebook": str(executed_path),
                "started": started,
                "finished": finished,
                "total_seconds": total,
                "peak_memory_gb": peak_gb,
                "summary": {"ok": n_ok, "error": n_err, "skipped": n_skip, "not_run": n_pending},
                "cells": results,
            },
            fh,
            indent=2,
            ensure_ascii=False,
        )

    lines = [
        f"# Notebook run report",
        "",
        f"- Notebook: `{args.notebook}`",
        f"- Executed copy (with all outputs): `{executed_path}`",
        f"- Started: {started}",
        f"- Finished: {finished or 'still running / interrupted'}",
        f"- Total cell runtime: {fmt_seconds(total)}",
        f"- Peak kernel memory: {f'{peak_gb:.1f} GB' if peak_gb else 'n/a'}",
        f"- Result: **{n_ok} ok, {n_err} errors, {n_skip} skipped, {n_pending} not run**",
        "",
        "Cell = position in the notebook counting from 1, including text cells (`--list` shows them).",
        "",
        "| Cell | Status | Time | Peak mem | First line | Error |",
        "|---:|---|---:|---:|---|---|",
    ]
    for r in results:
        err = (r["error"] or "").replace("|", "\\|").replace("\n", " ")[:200]
        fl = r["first_line"].replace("|", "\\|")
        mem = f"{r['peak_gb']:.1f} GB" if r.get("peak_gb") else ""
        lines.append(f"| {r['index']} | {r['status']} | {fmt_seconds(r['seconds'])} | {mem} | `{fl}` | {err} |")
    Path(report_md).write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--notebook", default="main_playground_clean.ipynb")
    p.add_argument("--run-dir", default=None,
                   help="folder for this run's report, executed notebook and log copies "
                        "(default: outputs/run_logs/<date_time>)")
    p.add_argument("--log-dir", default="outputs/logs",
                   help="folder whose logs are copied into the run folder at the end (only files written during this run)")
    p.add_argument("--skip", nargs="*", default=[], help="cell indices to skip, e.g. 2 44 49 or 40-45")
    p.add_argument("--only", nargs="*", default=None, help="run only these cell indices")
    p.add_argument("--continue-on-error", action="store_true", help="keep running after a cell fails")
    p.add_argument("--timeout", type=int, default=-1, help="per-cell timeout in seconds (-1 = no limit)")
    p.add_argument("--kernel", default="python3", help="Jupyter kernel name (default: python3 of this env)")
    p.add_argument("--list", action="store_true", help="list cells and exit")
    p.add_argument("--rerun-all", action="store_true",
                   help="reprocess every sample instead of loading already-processed samples from outputs/sample_cache")
    p.add_argument("--rerun", nargs="*", default=[], metavar="SAMPLE",
                   help="reprocess only these samples (8-character sample names, e.g. JMT49_CK)")
    args = p.parse_args()

    nb_path = Path(args.notebook)
    nb = nbformat.read(nb_path, as_version=4)

    if args.list:
        for i, c in enumerate(nb.cells):
            print(f"[{i + 1:3d}] {c.cell_type:8s} {first_line(c)}")
        return 0

    skip = parse_indices(args.skip)
    only = parse_indices(args.only) if args.only is not None else None

    # run with the working directory = folder containing the notebook (paths in the notebook are relative)
    # absolute() (not resolve()) keeps mapped drives / symlinks as given; resolve() can produce
    # paths longer than Windows' 260-character limit
    workdir = nb_path.absolute().parent
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_start = time.time()
    run_dir = Path(args.run_dir) if args.run_dir else Path("outputs/run_logs") / stamp
    if not run_dir.is_absolute():
        run_dir = workdir / run_dir
    run_dir.mkdir(parents=True, exist_ok=True)
    executed_path = run_dir / f"{nb_path.stem}_executed.ipynb"
    report_md = run_dir / "run_report.md"
    report_json = run_dir / "run_report.json"
    log_dir = Path(args.log_dir) if Path(args.log_dir).is_absolute() else workdir / args.log_dir
    # the notebook's first cell opens its log file in outputs/logs before it creates the folder tree
    log_dir.mkdir(parents=True, exist_ok=True)

    # headless plotting: plt.show() becomes a no-op, savefig still works
    os.environ["MPLBACKEND"] = "Agg"
    os.environ.setdefault("PYTHONIOENCODING", "utf-8")
    # sample cache controls (read by mayo/sample_cache.py inside the kernel)
    os.environ["MAYO_RERUN_ALL"] = "1" if args.rerun_all else ""
    os.environ["MAYO_RERUN_SAMPLES"] = " ".join(args.rerun)

    for c in nb.cells:
        if c.cell_type == "code":
            c.outputs = []
            c.execution_count = None

    results = [
        {"index": i + 1, "type": c.cell_type, "first_line": first_line(c), "status": "not run", "seconds": 0.0, "error": None, "peak_gb": None}
        for i, c in enumerate(nb.cells)
    ]
    results = [r for r in results if r["type"] == "code"]
    by_index = {r["index"]: r for r in results}

    started = datetime.now().isoformat(timespec="seconds")
    print(f"[{started}] Running {nb_path} in {workdir}")
    print(f"Executed notebook -> {executed_path}")
    print(f"Report            -> {report_md}")
    if args.rerun_all:
        print("Sample cache      -> ignored (--rerun-all): every sample is reprocessed")
    elif args.rerun:
        print(f"Sample cache      -> reprocessing {', '.join(args.rerun)}; other samples loaded from cache if present")
    else:
        print("Sample cache      -> only samples not processed before are run (outputs/sample_cache)")
    sys.stdout.flush()

    client = NotebookClient(
        nb,
        timeout=None if args.timeout < 0 else args.timeout,
        kernel_name=args.kernel,
        resources={"metadata": {"path": str(workdir)}},
        allow_errors=True,  # we inspect errors ourselves so we can decide to continue or stop
        record_timing=True,
    )

    exit_code = 0
    watcher = None
    try:
        with client.setup_kernel(cwd=str(workdir)):
            try:
                watcher = MemoryWatcher(client.km.provisioner.process.pid)
            except Exception:
                watcher = None
            for i, cell in enumerate(nb.cells):
                if cell.cell_type != "code":
                    continue
                n = i + 1  # cell number as counted in VS Code
                r = by_index[n]
                if n in skip or (only is not None and n not in only):
                    r["status"] = "skipped"
                    continue

                print(f"[{datetime.now():%H:%M:%S}] cell {n:3d} | {r['first_line'][:70]}", flush=True)
                if watcher:
                    watcher.start_cell()
                t0 = time.time()
                try:
                    client.execute_cell(cell, i)
                    err = error_text(cell)
                except (CellTimeoutError, DeadKernelError) as e:
                    err = f"{type(e).__name__}: {e}"
                except CellExecutionError as e:
                    err = str(e).splitlines()[-1] if str(e) else "CellExecutionError"
                r["seconds"] = time.time() - t0
                if watcher:
                    r["peak_gb"] = round(watcher.cell_peak, 2)

                if err:
                    r["status"] = "error"
                    r["error"] = err
                    print(f"           -> ERROR after {fmt_seconds(r['seconds'])}: {err}", flush=True)
                else:
                    r["status"] = "ok"
                    print(f"           -> ok ({fmt_seconds(r['seconds'])})", flush=True)

                nbformat.write(nb, executed_path)
                write_report(results, args, started, None, report_md, report_json, executed_path,
                             watcher.overall_peak if watcher else None)

                if err:
                    exit_code = 1
                    if isinstance(err, str) and err.startswith("DeadKernelError"):
                        print("Kernel died (often out of memory). Stopping.", flush=True)
                        break
                    if not args.continue_on_error:
                        print("Stopping at first error (use --continue-on-error to keep going).", flush=True)
                        break
    except Exception:
        traceback.print_exc()
        exit_code = 2
    finally:
        finished = datetime.now().isoformat(timespec="seconds")
        if watcher:
            watcher.stop()
        nbformat.write(nb, executed_path)
        write_report(results, args, started, finished, report_md, report_json, executed_path,
                     watcher.overall_peak if watcher else None)
        copied = archive_logs(log_dir, run_dir / "logs", since=run_start)
        if copied:
            print(f"Copied {len(copied)} log file(s) written during this run to {run_dir / 'logs'}", flush=True)

    n_ok = sum(r["status"] == "ok" for r in results)
    n_err = sum(r["status"] == "error" for r in results)
    peak = f", peak memory {watcher.overall_peak:.1f} GB" if watcher and watcher.overall_peak else ""
    print(f"[{finished}] Done: {n_ok} ok, {n_err} errors{peak}. Report: {report_md}", flush=True)
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
