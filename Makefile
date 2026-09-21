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
	@echo "  make review-sex-qc Review sex-QC results and continue"
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
	candidate="results/review/sex/candidate_sex_thresholds.yaml"; \
	accepted="config/accepted_sex_thresholds.yaml"; \
	if [[ ! -s "$$summary" ]]; then \
		echo "ERROR: sex-QC summary not found:" >&2; \
		echo "  $$summary" >&2; \
		echo "Run 'make run' first." >&2; \
		exit 1; \
	fi; \
	if [[ ! -s "$$candidate" ]]; then \
		echo "ERROR: candidate thresholds not found:" >&2; \
		echo "  $$candidate" >&2; \
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
	echo "Candidate thresholds:"; \
	echo "----------------------------------------"; \
	cat "$$candidate"; \
	echo ""; \
	printf "Continue with these thresholds? [y/n]: "; \
	read -r answer; \
	case "$$answer" in \
		y|Y|yes|YES) \
			cp "$$candidate" "$$accepted"; \
			echo ""; \
			echo "Thresholds accepted:"; \
			cat "$$accepted"; \
			echo ""; \
			echo "Continuing the workflow..."; \
			$(MAKE) run \
			;; \
		n|N|no|NO) \
			echo ""; \
			printf "New maximum F for females: "; \
			read -r female_max; \
			printf "New minimum F for males: "; \
			read -r male_min; \
			if ! awk \
				-v female="$$female_max" \
				-v male="$$male_min" \
				'BEGIN { \
					number = "^-?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$$"; \
					if (female !~ number || male !~ number) exit 1; \
					if ((female + 0) >= (male + 0)) exit 1; \
				}'; then \
				echo "ERROR: thresholds must be numeric and" >&2; \
				echo "female_max_f must be lower than male_min_f." >&2; \
				exit 1; \
			fi; \
			printf "female_max_f: %s\nmale_min_f: %s\n" \
				"$$female_max" \
				"$$male_min" \
				> "$$accepted"; \
			echo ""; \
			echo "New thresholds saved:"; \
			cat "$$accepted"; \
			echo ""; \
			echo "Continuing the workflow..."; \
			$(MAKE) run \
			;; \
		*) \
			echo ""; \
			echo "No valid choice entered."; \
			echo "The workflow was not continued."; \
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