"""Shared configuration values and helper functions."""

from pathlib import Path


_MISSING = object()


def cfg(key_path, default=_MISSING):
    """Return a config value addressed as ``section/key``."""
    value = config
    for key in key_path.split("/"):
        if not isinstance(value, dict) or key not in value:
            if default is not _MISSING:
                return default
            raise WorkflowError(
                f"Missing required configuration entry: '{key_path}'"
            )
        value = value[key]
    return value


RESULTS_DIR = Path(cfg("output_dir", "results"))


def ws_path(relative_path):
    """Return a path below the configured workflow output directory."""
    path = Path(relative_path)
    if path.is_absolute():
        raise WorkflowError(
            "ws_path() requires a relative path; received "
            f"'{relative_path}'"
        )
    return str(RESULTS_DIR / path)
