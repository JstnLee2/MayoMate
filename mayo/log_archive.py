"""Keep a dated copy of a log file so earlier runs are not lost when the log is overwritten."""
import os
import shutil
from datetime import datetime

ARCHIVE_DIR = "outputs/run_logs/scripts"


def keep_log_copy(log_path: str, archive_dir: str = ARCHIVE_DIR) -> str | None:
    """
    Copy log_path to <archive_dir>/<log name>_<YYYYmmdd_HHMMSS>.log and return the new path.
    The log in outputs/logs/ always holds the latest run; the dated copies keep every run.
    """
    import logging
    if not os.path.isfile(log_path):
        return None
    for handler in logging.getLogger().handlers:  # make sure everything is written before copying
        try:
            handler.flush()
        except Exception:
            pass
    os.makedirs(archive_dir, exist_ok=True)
    stem, ext = os.path.splitext(os.path.basename(log_path))
    dest = os.path.join(archive_dir, f"{stem}_{datetime.now():%Y%m%d_%H%M%S}{ext or '.log'}")
    shutil.copy2(log_path, dest)
    return dest
