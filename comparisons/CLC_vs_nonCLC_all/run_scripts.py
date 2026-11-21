"""Run cluster_analysis.py, cluster_coordination.py and cluster_association_chi2.py in three ways
after the notebook has finished, and keep each run's outputs separately in script_runs/<run>/.

  per genotype       script_runs/pair_<genotype>/   <genotype>CLC (reference) vs <genotype> (own caller)
  CLC set only       script_runs/CLC_set/           every CLC genotype vs ung1∆CLC
  own-caller set     script_runs/own_set/           every own-caller genotype vs ung1∆

The scripts overwrite some of their outputs (figures, cluster_length_meta_df_JT.csv) on every run,
so every file in outputs/ and figures/ written during a run is copied into that run's folder.
Usage:  ..\\..\\myenv\\Scripts\\python.exe run_scripts.py
"""
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

PY = sys.executable
SCRIPTS = ["cluster_analysis.py", "cluster_coordination.py", "cluster_association_chi2.py"]
GENOTYPES = sorted({p["genotype"] for p in json.load(open("pairs.json", encoding="utf-8"))["pairs"].values()})
REF = "ung1∆"

runs = [(f"pair_{g}", g + "CLC", [g + "CLC", g]) for g in GENOTYPES]
runs += [("CLC_set", REF + "CLC", [g + "CLC" for g in GENOTYPES]),
         ("own_set", REF, list(GENOTYPES))]

env = dict(os.environ, MPLBACKEND="Agg", PYTHONIOENCODING="utf-8")
summary = []
for name, reference, genotypes in runs:
    safe = re.sub(r"[^\w∆.-]+", "_", name)
    dest = Path("script_runs") / safe
    dest.mkdir(parents=True, exist_ok=True)
    for script in SCRIPTS:
        t0 = time.time()
        cmd = [PY, script, "--reference", reference, "--genotypes", *genotypes]
        with open(dest / f"{script[:-3]}_stdout.txt", "w", encoding="utf-8") as fh:
            rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env).returncode
        copied = 0
        for top in ["outputs", "figures"]:
            for root, _, files in os.walk(top):
                if "run_logs" in root or "sample_cache" in root or "partial_pickles" in root:
                    continue
                for f in files:
                    src = Path(root) / f
                    if src.stat().st_mtime >= t0:
                        # \\?\ prefix: some figure names push the copy past Windows' 260-character limit
                        out = Path("\\\\?\\" + os.path.abspath(dest / src))
                        out.parent.mkdir(parents=True, exist_ok=True)
                        shutil.copy2("\\\\?\\" + os.path.abspath(src), out)
                        copied += 1
        line = f"{name:40s} {script:28s} exit {rc}  {time.time() - t0:6.0f}s  {copied} files"
        print(line, flush=True)
        summary.append(dict(run=name, script=script, reference=reference, genotypes=genotypes, exit_code=rc, files=copied))
json.dump(summary, open("script_runs/summary.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
