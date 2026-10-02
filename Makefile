SHELL := /bin/bash

SNAKEMAKE ?= snakemake
SNAKEFILE ?= Snakefile
CONFIGFILE ?= config/config.yaml
PROFILE ?= profiles/slurm
CORES ?= 1

COMMON_ARGS = --snakefile $(SNAKEFILE) --configfile $(CONFIGFILE)

.DEFAULT_GOAL := help

.PHONY: \
	help \
	dry-run \
	run \
	run-local \
	review-sex-qc \
	lint \
	summary \
	dag \
	rulegraph \
	unlock


help:
	@echo "Available targets:"
	@echo "  make dry-run       Validate the DAG without running jobs"
	@echo "  make run           Submit jobs with the SLURM profile"
	@echo "  make run-local     Run locally (override with CORES=N)"
	@echo "  make review-sex-qc Review sex-QC thresholds after a completed run"
	@echo "  make lint          Run the Snakemake workflow linter"
	@echo "  make summary       Show the status of output files"
	@echo "  make dag           Write the job DAG to dag.svg"
	@echo "  make rulegraph     Write the rule graph to rulegraph.svg"
	@echo "  make unlock        Remove a stale Snakemake lock"


dry-run:
	$(SNAKEMAKE) $(COMMON_ARGS) \
		--dry-run \
		--printshellcmds \
		--cores 1


run:
	$(SNAKEMAKE) $(COMMON_ARGS) --profile $(PROFILE)


run-local:
	$(SNAKEMAKE) $(COMMON_ARGS) \
		--cores $(CORES) \
		--software-deployment-method conda \
		--printshellcmds


review-sex-qc:
	@set -euo pipefail; \
	summary="results/review/sex/sex_qc.summary.tsv"; \
	plot="results/review/sex/sex_F_distribution.pdf"; \
	table="results/review/sex/sex_candidate_classification.tsv"; \
	if [[ ! -s "$$summary" ]]; then \
		echo "ERROR: sex-QC summary not found:" >&2; \
		echo "  $$summary" >&2; \
		echo "Run 'make run' first." >&2; \
		exit 1; \
	fi; \
	echo ""; \
	echo "========================================"; \
	echo "Sex-QC manual review"; \
	echo "========================================"; \
	echo ""; \
	echo "Files to inspect:"; \
	echo "  Plot: $$plot"; \
	echo "  Table: $$table"; \
	echo ""; \
	echo "Summary:"; \
	echo "----------------------------------------"; \
	if command -v column >/dev/null 2>&1; then \
		column -t -s $$'\t' "$$summary"; \
	else \
		cat "$$summary"; \
	fi; \
	echo ""; \
	printf "Accept the sex-QC thresholds used in this run? [y/n]: "; \
	read -r answer; \
	case "$$answer" in \
		y|Y|yes|YES) \
			echo ""; \
			echo "Sex-QC thresholds accepted."; \
			echo "Existing pipeline results are kept unchanged."; \
			;; \
		n|N|no|NO) \
			echo ""; \
			printf "New maximum F for females: "; \
			read -r female_max; \
			printf "New minimum F for males: "; \
			read -r male_min; \
			python3 workflow/scripts/update_sex_thresholds.py \
				--config "$(CONFIGFILE)" \
				--female-max-f "$$female_max" \
				--male-min-f "$$male_min"; \
			echo ""; \
			echo "Thresholds changed. Rerunning sex-QC-dependent workflow steps..."; \
			$(MAKE) run; \
			;; \
		*) \
			echo ""; \
			echo "ERROR: enter y or n." >&2; \
			exit 1 \
			;; \
	esac

lint:
	$(SNAKEMAKE) $(COMMON_ARGS) --lint


summary:
	$(SNAKEMAKE) $(COMMON_ARGS) --summary


dag:
	$(SNAKEMAKE) $(COMMON_ARGS) --dag | dot -Tsvg > dag.svg


rulegraph:
	$(SNAKEMAKE) $(COMMON_ARGS) \
		--rulegraph | dot -Tsvg > rulegraph.svg


unlock:
	$(SNAKEMAKE) $(COMMON_ARGS) --unlock