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
