"""
Per-sample cache: re-running the notebook only processes samples that are new (or whose input
file changed). Results for samples that were already processed are loaded from
outputs/sample_cache/<stage>/.

Only the slow per-sample steps are cached (import + filters + ploidy, recombination maps,
JT clusters). Steps that compare samples with each other (e.g. filterCommonMutations) and all
downstream analyses/figures still run on the full set of samples every time.

A sample is processed again when:
  - its input file is new, or its size / modification time changed,
  - the parameters passed for that stage changed (e.g. SNP_NUM, filter thresholds),
  - the parental reference (parental_origin_db_df) changed,
  - you force it: environment variable MAYO_RERUN_ALL=1 (everything) or
    MAYO_RERUN_SAMPLES="JMT49_CK JMT50_CK" (only those). run_notebook.py sets these with
    --rerun-all / --rerun.
Changes to the mayo code itself are NOT detected; use --rerun-all after changing how samples are
processed (or delete outputs/sample_cache).
"""

import hashlib
import json
import logging
import os
import pickle
import re
from pathlib import Path
from typing import Callable

import pandas as pd

CACHE_ROOT = Path("outputs/sample_cache")
CACHE_VERSION = 1  # bump to invalidate every cache entry


def rerun_all() -> bool:
    return os.environ.get("MAYO_RERUN_ALL", "").strip().lower() in ("1", "true", "yes")


def rerun_samples() -> set:
    return set(re.split(r"[,\s]+", os.environ.get("MAYO_RERUN_SAMPLES", "").strip())) - {""}


def _forced(name: str) -> bool:
    return rerun_all() or name in rerun_samples()


def make_key(*parts) -> str:
    return hashlib.sha1(json.dumps(parts, sort_keys=True, default=str).encode("utf-8")).hexdigest()


def file_fingerprint(path) -> str:
    st = os.stat(path)
    return f"{os.path.basename(path)}|{st.st_size}|{st.st_mtime_ns}"


def dir_fingerprint(directory, ext: str | None = None) -> list:
    return [file_fingerprint(os.path.join(directory, f)) for f in sorted(os.listdir(directory))
            if os.path.isfile(os.path.join(directory, f)) and (ext is None or f.endswith(ext))]


def df_fingerprint(df: pd.DataFrame) -> str:
    return hashlib.sha1(pd.util.hash_pandas_object(df, index=False).values.tobytes()).hexdigest()


class StageCache:
    """One pickle per sample per stage: {"key": <fingerprint>, "obj": <payload>}."""

    def __init__(self, stage: str):
        self.stage = stage
        self.dir = CACHE_ROOT / stage

    def _path(self, entry_id: str) -> Path:
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", entry_id)
        if safe != entry_id or len(safe) > 60:
            safe = safe[:40] + "_" + hashlib.sha1(entry_id.encode()).hexdigest()[:10]
        return self.dir / f"{safe}.pkl"

    def get(self, entry_id: str, key: str, name: str | None = None):
        if _forced(name or entry_id):
            return None
        path = self._path(entry_id)
        if not path.exists():
            return None
        try:
            with open(path, "rb") as fh:
                payload = pickle.load(fh)
        except Exception as e:
            logging.warning(f"[cache:{self.stage}] could not read {path} ({e}); recomputing")
            return None
        return payload["obj"] if payload.get("key") == key else None

    def put(self, entry_id: str, key: str, obj) -> None:
        self.dir.mkdir(parents=True, exist_ok=True)
        path = self._path(entry_id)
        tmp = path.with_suffix(".tmp")
        with open(tmp, "wb") as fh:
            pickle.dump({"key": key, "obj": obj}, fh, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp, path)


def _report(stage: str, n_cached: int, new_names: list) -> None:
    msg = f"[cache:{stage}] {n_cached} sample(s) loaded from cache, {len(new_names)} processed"
    if new_names:
        shown = ", ".join(new_names[:20]) + (" ..." if len(new_names) > 20 else "")
        msg += f": {shown}"
    print(msg)
    logging.info(msg)


def createParentalOriginDF_cached(parents_path: str, **kwargs) -> pd.DataFrame:
    """createParentalOriginDF(), skipped when the parent files and settings are unchanged."""
    from .base_functions import createParentalOriginDF
    from .settings import config as cfg

    settings_dir = Path(__file__).parent / "settings"
    key = make_key(CACHE_VERSION, dir_fingerprint(parents_path), kwargs,
                   file_fingerprint(settings_dir / "problematic_positions.py"),
                   file_fingerprint(cfg["repetitive_regions_fasta_path"]))
    cache = StageCache("parental_origin")
    df = cache.get("parental_origin_db_df", key, name="__parents__")
    save_path = kwargs.get("parental_origin_db_df_path_save")
    if df is not None:
        print("[cache:parental_origin] parental reference unchanged, loaded from cache")
        if save_path and not os.path.exists(save_path):
            df.to_csv(save_path, sep="\t", index=False)
        return df
    df = createParentalOriginDF(parents_path, **kwargs)
    cache.put("parental_origin_db_df", key, df)
    return df


def import_samples_cached(
        directory: str,
        kind: str,
        stage: str,
        process_new: Callable[[list], None] | None = None,
        params: dict | None = None,
        genotype_dict: dict | None = None,
        format: str = "csv",
        log: bool = False) -> tuple[list, list]:
    """csvDataFrameImport() + per-sample processing, but only for samples not yet in the cache.

    process_new(list_of_new_samples) must modify the samples in place; the processed samples are
    then saved to the cache. Returns (all_samples_sorted, names_of_newly_processed_samples).
    Every sample gets a `cache_key` attribute that later cached stages build on.
    """
    from .base_functions import read_csv_file, create_table_object

    cache = StageCache(stage)
    samples, todo = [], []
    for filename in sorted(os.listdir(directory)):
        if not filename.endswith(f".{format}"):
            continue
        name = filename[:8]
        if genotype_dict is not None and name not in genotype_dict:
            logging.warning(f"Sample {name} not found in genotype_dict. Skipping...")
            continue
        path = os.path.join(directory, filename)
        key = make_key(CACHE_VERSION, stage, file_fingerprint(path), params)
        obj = cache.get(filename, key, name=name)
        if obj is not None:
            obj.cache_key = key
            samples.append(obj)
        else:
            todo.append((filename, path, name, key))

    new = []
    for filename, path, name, key in todo:
        table = create_table_object(read_csv_file(path, format), kind, name)
        if log:
            table.printRows()
        table.extractSNVs(log=log)
        table.cache_key = key
        table.cache_file = filename
        new.append(table)

    if new:
        if process_new is not None:
            process_new(new)
        for table in new:
            cache.put(table.cache_file, table.cache_key, table)

    _report(stage, len(samples), [t.name for t in new])
    samples.extend(new)
    samples.sort()
    logging.info(f"found {len(samples)} samples in {directory}")
    return samples, [t.name for t in new]


def createRecombinationMaps_cached(
        samples_parental_SNPs: list,
        parental_origin_db_df: pd.DataFrame,
        df_chr_lengths: pd.DataFrame,
        starts,
        draw_map: bool = False,
        saveMap: bool = True,
        min_SNPs_switch_support: int = 1,
        snps_and_rec_ev_file: str = "outputs/parent_snps_and_rec_ev.txt") -> list:
    """createRecombinationMaps(), run only for samples not yet in the cache.

    Like the original, it updates the sample objects: samples_parental_SNPs is changed in place
    (cached samples are replaced by their cached post-recombination-map version).
    Returns the switches list (same order as samples_parental_SNPs).
    """
    from .base_functions import createRecombinationMaps

    stage = f"recombination_{min_SNPs_switch_support}snp"
    cache = StageCache(stage)
    ref_fp = df_fingerprint(parental_origin_db_df)
    results, todo = {}, []
    for idx, s in enumerate(samples_parental_SNPs):
        entry_id = getattr(s, "cache_file", s.name)
        key = make_key(CACHE_VERSION, stage, getattr(s, "cache_key", None), ref_fp)
        hit = cache.get(entry_id, key, name=s.name) if getattr(s, "cache_key", None) else None
        if hit is not None:
            results[idx] = hit
        else:
            todo.append((idx, s, entry_id, key))

    if todo:
        tmp_file = snps_and_rec_ev_file + ".new_samples.tmp"
        new_switches = createRecombinationMaps(
            [s for _, s, _, _ in todo], parental_origin_db_df, df_chr_lengths, starts,
            draw_map=draw_map, saveMap=saveMap, min_SNPs_switch_support=min_SNPs_switch_support,
            snps_and_rec_ev_file=tmp_file)
        if os.path.exists(tmp_file):
            os.remove(tmp_file)
        for (idx, s, entry_id, key), sw in zip(todo, new_switches):
            results[idx] = (s, sw)
            if getattr(s, "cache_key", None):
                cache.put(entry_id, key, (s, sw))

    _report(stage, len(samples_parental_SNPs) - len(todo), [s.name for _, s, _, _ in todo])
    ordered = [results[i] for i in range(len(samples_parental_SNPs))]
    samples_parental_SNPs[:] = [s for s, _ in ordered]

    # the combined text file is rewritten for ALL samples
    with open(snps_and_rec_ev_file, "w") as f:
        for s in samples_parental_SNPs:
            f.write(getattr(s, "output_str", ""))

    return [sw for _, sw in ordered]
