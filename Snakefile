configfile: "config/config.yaml"

include: "workflow/rules/common.smk"
include: "workflow/rules/standardization.smk"
include: "workflow/rules/sample_qc.smk"
include: "workflow/rules/structure.smk"
include: "workflow/rules/final.smk"

