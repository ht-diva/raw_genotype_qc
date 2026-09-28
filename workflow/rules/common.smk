"""Shared configuration values and helper functions."""

from pathlib import Path


# Accept the sample_metadata section from existing BELIEVE configurations.
if "phenotype" not in config and "sample_metadata" in config:
    config["phenotype"] = dict(config["sample_metadata"])
if isinstance(config.get("phenotype"), dict):
    for target, legacy in (("male_values", "male_values"),
                           ("female_values", "female_value")):
        if target not in config["phenotype"] and legacy in config:
            config["phenotype"][target] = config[legacy]

# A second inclusion mechanism could unexpectedly restrict the disabled mode.
if config.get("external_samples", {}).get("keep_files"):
    raise WorkflowError(
        "Move external_samples/keep_files to keep_ids/proteomics or "
        "keep_ids/metabolomics; keep_ids/enabled controls initial inclusion."
    )


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


# Initial sample selection: keep IDs from proteomics or metabolomics -----------

def selected_keep_ids_files():
    enabled = cfg("keep_ids/enabled")

    if not isinstance(enabled, bool):
        raise WorkflowError("keep_ids/enabled must be true or false")
    if not enabled:
        return []

    assay = cfg("keep_ids/type")
    if assay not in ("proteomics", "metabolomics"):
        raise WorkflowError(
            "keep_ids/type must be 'proteomics' or 'metabolomics'"
        )
    path = cfg(f"keep_ids/{assay}")
    if not isinstance(path, str) or not path.strip():
        raise WorkflowError(
            f"keep_ids/{assay} must contain an ID file path when enabled"
        )
    return [path]


def external_remove_files():
    """External exclusions apply independently of the omics keep switch."""
    files = cfg("external_samples/remove_files", [])
    if not isinstance(files, list):
        raise WorkflowError("external_samples/remove_files must be a list")
    path = cfg("external_samples_exclusions", None)
    if path is not None:
        if not isinstance(path, str) or not path.strip():
            raise WorkflowError("external_samples_exclusions must be a path or null")
        files = [path] + files
    if any(not isinstance(p, str) or not p.strip() for p in files):
        raise WorkflowError("External removal files must contain nonempty paths")
    return list(dict.fromkeys(files))


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
