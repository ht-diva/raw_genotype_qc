#!/usr/bin/env bash

# Inspect a merged PLINK BED dataset before downstream quality control.
#
# This script:
#   1. Verifies that the BED, BIM, and FAM input files exist and are not empty.
#   2. Counts the number of samples and variants.
#   3. Checks that BIM and FAM rows contain at least six fields.
#   4. Identifies duplicated FID/IID sample identifiers.
#   5. Identifies duplicated variant IDs.
#   6. Identifies variants sharing the same chromosome and position.
#   7. Creates a chromosome-name mapping for PLINK --rename-chrs.
#   8. Counts missing variant IDs and summarizes reported sample sex.
#   9. Writes all results to a summary TSV file.
#
# The script does not modify the original BED, BIM, or FAM files.

set -euo pipefail


# ----------------------------------------------------------
# 1. Initialize command-line arguments
# ----------------------------------------------------------

bed=""
bim=""
fam=""
chrmap=""
summary=""
dup_samples=""
dup_ids=""
dup_pos=""
fail_dup="True"


# ----------------------------------------------------------
# 2. Parse command-line arguments
# ----------------------------------------------------------

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bed)
            bed="$2"
            shift 2
            ;;
        --bim)
            bim="$2"
            shift 2
            ;;
        --fam)
            fam="$2"
            shift 2
            ;;
        --chromosome-map)
            chrmap="$2"
            shift 2
            ;;
        --summary)
            summary="$2"
            shift 2
            ;;
        --duplicate-samples)
            dup_samples="$2"
            shift 2
            ;;
        --duplicate-variant-ids)
            dup_ids="$2"
            shift 2
            ;;
        --duplicate-positions)
            dup_pos="$2"
            shift 2
            ;;
        --fail-duplicate-samples)
            fail_dup="$2"
            shift 2
            ;;
        *)
            echo "ERROR: unknown argument: $1" >&2
            exit 2
            ;;
    esac
done


# ----------------------------------------------------------
# 3. Verify that the PLINK input files exist and are not empty
# ----------------------------------------------------------

for input_file in "$bed" "$bim" "$fam"; do
    if [[ ! -s "$input_file" ]]; then
        echo \
            "ERROR: missing or empty input file: $input_file" \
            >&2
        exit 1
    fi
done

# Create the output directory if it does not already exist.
mkdir -p "$(dirname "$summary")"


# ----------------------------------------------------------
# 4. Count samples and variants
# ----------------------------------------------------------

# Each non-empty FAM row represents one sample.
n_samples=$(
    awk '
        NF {
            n++
        }
        END {
            print n + 0
        }
    ' "$fam"
)

# Each non-empty BIM row represents one variant.
n_variants=$(
    awk '
        NF {
            n++
        }
        END {
            print n + 0
        }
    ' "$bim"
)

# Stop if the dataset contains no samples or no variants.
if [[ "$n_samples" -eq 0 || "$n_variants" -eq 0 ]]; then
    echo \
        "ERROR: BED fileset contains zero samples or variants" \
        >&2
    exit 1
fi


# ----------------------------------------------------------
# 5. Check the minimum number of BIM and FAM fields
# ----------------------------------------------------------

# A valid FAM row must contain at least:
# FID IID PAT MAT SEX PHENOTYPE
awk '
    NF && NF < 6 {
        bad = 1
    }
    END {
        exit bad
    }
' "$fam" || {
    echo "ERROR: FAM row has fewer than 6 fields" >&2
    exit 1
}

# A valid BIM row must contain at least:
# CHR ID CM POS A1 A2
awk '
    NF && NF < 6 {
        bad = 1
    }
    END {
        exit bad
    }
' "$bim" || {
    echo "ERROR: BIM row has fewer than 6 fields" >&2
    exit 1
}


# ----------------------------------------------------------
# 6. Identify duplicated sample identifiers
# ----------------------------------------------------------

# Samples are considered duplicated when the same FID/IID
# combination occurs more than once in the FAM file.
{
    printf 'FID\tIID\tN\n'

    awk '
        NF {
            print $1 "\t" $2
        }
    ' "$fam" |
        sort |
        uniq -c |
        awk '
            $1 > 1 {
                print $2 "\t" $3 "\t" $1
            }
        '
} > "$dup_samples"

n_dup_samples=$(
    awk '
        NR > 1 {
            n++
        }
        END {
            print n + 0
        }
    ' "$dup_samples"
)

# Optionally stop the workflow when duplicated samples are found.
# Accepted true values are: true, 1, yes, and y.
if (
    [[ "$n_dup_samples" -gt 0 ]] &&
    [[ "${fail_dup,,}" =~ ^(true|1|yes|y)$ ]]
); then
    echo \
        "ERROR: $n_dup_samples duplicate FID/IID sample IDs found" \
        >&2
    exit 1
fi


# ----------------------------------------------------------
# 7. Identify duplicated variant IDs
# ----------------------------------------------------------

# The variant ID is stored in the second column of the BIM file.
{
    printf 'ID\tN\n'

    awk '
        NF {
            print $2
        }
    ' "$bim" |
        sort |
        uniq -c |
        awk '
            $1 > 1 {
                print $2 "\t" $1
            }
        '
} > "$dup_ids"

n_dup_ids=$(
    awk '
        NR > 1 {
            n++
        }
        END {
            print n + 0
        }
    ' "$dup_ids"
)


# ----------------------------------------------------------
# 8. Identify duplicated genomic positions
# ----------------------------------------------------------

# Chromosome labels are normalized before checking positions:
#   chr1, CHR1, Chr1 -> 1
#   23               -> X
#   24               -> Y
#   25               -> XY
#   26 or M           -> MT
#
# Variants are considered to share a position when they have
# the same normalized chromosome and base-pair coordinate.
awk '
    function normalize_chromosome(chromosome, normalized) {
        sub(/^chr/, "", chromosome)
        sub(/^CHR/, "", chromosome)
        sub(/^Chr/, "", chromosome)

        normalized = toupper(chromosome)

        if (normalized == "23") {
            normalized = "X"
        } else if (normalized == "24") {
            normalized = "Y"
        } else if (normalized == "25") {
            normalized = "XY"
        } else if (normalized == "26" || normalized == "M") {
            normalized = "MT"
        }

        return normalized
    }

    NF {
        print normalize_chromosome($1) "\t" $4
    }
' "$bim" |
    sort -k1,1V -k2,2n |
    uniq -c |
    awk '
        BEGIN {
            print "CHR\tPOS\tN"
        }

        $1 > 1 {
            print $2 "\t" $3 "\t" $1
        }
    ' > "$dup_pos"

n_dup_pos=$(
    awk '
        NR > 1 {
            n++
        }
        END {
            print n + 0
        }
    ' "$dup_pos"
)


# ----------------------------------------------------------
# 9. Create a chromosome-renaming map
# ----------------------------------------------------------

# The output contains one row for every unique original
# chromosome label:
#
# original_label    standardized_label
#
# The file has no header so that it can be used directly with:
# plink --rename-chrs chromosome_map.tsv
awk '
    function normalize_chromosome(chromosome, normalized) {
        sub(/^chr/, "", chromosome)
        sub(/^CHR/, "", chromosome)
        sub(/^Chr/, "", chromosome)

        normalized = toupper(chromosome)

        if (normalized == "23") {
            normalized = "X"
        } else if (normalized == "24") {
            normalized = "Y"
        } else if (normalized == "25") {
            normalized = "XY"
        } else if (normalized == "26" || normalized == "M") {
            normalized = "MT"
        }

        return normalized
    }

    NF && !seen[$1]++ {
        print $1 "\t" normalize_chromosome($1)
    }
' "$bim" > "$chrmap"


# ----------------------------------------------------------
# 10. Summarize chromosome labels
# ----------------------------------------------------------

# Create a comma-separated list of the original chromosome labels.
raw_chrs=$(
    awk '
        !seen[$1]++ {
            if (chromosomes == "") {
                chromosomes = $1
            } else {
                chromosomes = chromosomes "," $1
            }
        }

        END {
            print chromosomes
        }
    ' "$bim"
)

# Create a comma-separated list of standardized chromosome labels.
std_chrs=$(
    awk '
        !seen[$2]++ {
            if (chromosomes == "") {
                chromosomes = $2
            } else {
                chromosomes = chromosomes "," $2
            }
        }

        END {
            print chromosomes
        }
    ' "$chrmap"
)


# ----------------------------------------------------------
# 11. Count missing variant IDs
# ----------------------------------------------------------

# Variant IDs represented by "." or "0" are treated as missing.
missing_ids=$(
    awk '
        NF && ($2 == "." || $2 == "0") {
            n++
        }

        END {
            print n + 0
        }
    ' "$bim"
)


# ----------------------------------------------------------
# 12. Summarize reported sample sex
# ----------------------------------------------------------

# PLINK FAM sex codes:
#   1 = male
#   2 = female
#   other values = unknown
male=$(
    awk '
        NF && $5 == 1 {
            n++
        }

        END {
            print n + 0
        }
    ' "$fam"
)

female=$(
    awk '
        NF && $5 == 2 {
            n++
        }

        END {
            print n + 0
        }
    ' "$fam"
)

unknown=$((n_samples - male - female))


# ----------------------------------------------------------
# 13. Write the final input summary
# ----------------------------------------------------------

{
    printf 'metric\tvalue\n'
    printf 'n_samples\t%s\n' "$n_samples"
    printf 'n_variants\t%s\n' "$n_variants"
    printf 'n_duplicate_sample_ids\t%s\n' "$n_dup_samples"
    printf 'n_duplicate_variant_ids\t%s\n' "$n_dup_ids"
    printf 'n_duplicate_chr_pos\t%s\n' "$n_dup_pos"
    printf 'n_missing_variant_ids\t%s\n' "$missing_ids"
    printf 'fam_sex_male\t%s\n' "$male"
    printf 'fam_sex_female\t%s\n' "$female"
    printf 'fam_sex_unknown\t%s\n' "$unknown"
    printf 'chromosomes_raw\t%s\n' "$raw_chrs"
    printf 'chromosomes_standardized\t%s\n' "$std_chrs"
} > "$summary"
