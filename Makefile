SNAKEMAKE ?= snakemake
SNAKEFILE ?= Snakefile
CONFIGFILE ?= config/config.yaml
PROFILE ?= profiles/slurm
CORES ?= 1

COMMON_ARGS = --snakefile $(SNAKEFILE) --configfile $(CONFIGFILE)

.DEFAULT_GOAL := help

.PHONY: help dry-run run run-local lint summary dag rulegraph unlock

help:
	@echo "Available targets:"
	@echo "  make dry-run    Validate the DAG without running jobs"
	@echo "  make run        Submit jobs with the SLURM profile"
	@echo "  make run-local  Run locally (override with CORES=N)"
	@echo "  make lint       Run the Snakemake workflow linter"
	@echo "  make summary    Show the status of output files"
	@echo "  make dag        Write the job DAG to dag.svg (Graphviz required)"
	@echo "  make rulegraph  Write the rule graph to rulegraph.svg (Graphviz required)"
	@echo "  make unlock     Remove a stale Snakemake working-directory lock"

dry-run:
	$(SNAKEMAKE) $(COMMON_ARGS) --dry-run --printshellcmds --cores 1

run:
	$(SNAKEMAKE) $(COMMON_ARGS) --profile $(PROFILE)

run-local:
	$(SNAKEMAKE) $(COMMON_ARGS) --cores $(CORES) \
		--software-deployment-method conda --printshellcmds

lint:
	$(SNAKEMAKE) $(COMMON_ARGS) --lint

summary:
	$(SNAKEMAKE) $(COMMON_ARGS) --summary

dag:
	$(SNAKEMAKE) $(COMMON_ARGS) --dag | dot -Tsvg > dag.svg

rulegraph:
	$(SNAKEMAKE) $(COMMON_ARGS) --rulegraph | dot -Tsvg > rulegraph.svg

unlock:
	$(SNAKEMAKE) $(COMMON_ARGS) --unlock
