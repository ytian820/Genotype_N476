version 1.0

# Purpose: filter samples by call rate, then estimate KING-robust kinship.
# Phenotype: none. Unit: one cohort/sample_set.
# Inputs: autosomal PGEN/PVAR/PSAM from 04; hardcalls, chromosomes 1-22 only.
# Outputs: sample-filtered trio, KING table and two PLINK logs.
# Retain call rate >= 98% (--mind 0.02); preserve input FID/IID and all variants.
# No MAF, HWE, LD pruning or automatic removal of relatives.

# Part 1. Workflow interface
workflow SampleCallRateAndKing {
    input {
        File pgen
        File pvar
        File psam
    }

    call FilterSamplesAndRunKing {
        input: pgen = pgen, pvar = pvar, psam = psam
    }

    output {
        File filtered_pgen = FilterSamplesAndRunKing.filtered_pgen
        File filtered_pvar = FilterSamplesAndRunKing.filtered_pvar
        File filtered_psam = FilterSamplesAndRunKing.filtered_psam
        File king_table = FilterSamplesAndRunKing.king_table
        Array[File] plink_logs = FilterSamplesAndRunKing.plink_logs
    }
}

# Part 2. Sequential sample filtering and relatedness estimation
task FilterSamplesAndRunKing {
    input {
        File pgen
        File pvar
        File psam
    }

    command <<<
        set -euo pipefail
        echo "Purpose: sample call-rate filtering followed by KING-robust"
        echo "Threshold: sample hardcall missingness <= 0.02 (call rate >= 98%)"
        echo "Variants: all autosomal sites from 04; no MAF filter or LD pruning"
        echo "Relatedness: report all sample pairs; do not remove relatives"
        echo "Resources: 1 CPU, 4 GB memory, 20 GB local disk"

        # Preflight: require PLINK2 and all three nonempty input files.
        command -v plink2 >/dev/null
        for input_file in "~{pgen}" "~{pvar}" "~{psam}"; do
            if [[ ! -s "$input_file" ]]; then
                echo "ERROR: an input PGEN/PVAR/PSAM file is missing or empty" >&2
                exit 1
            fi
        done

        ## Step 1. Remove only samples with missingness > 2%.
        plink2 --pgen "~{pgen}" --pvar "~{pvar}" --psam "~{psam}" \
            --mind 0.02 --make-pgen \
            --threads 1 --memory 3000 \
            --out autosomes.callrate98

        ## Step 2. Calculate KING-robust for all remaining sample pairs.
        plink2 --pfile autosomes.callrate98 \
            --make-king-table \
            --threads 1 --memory 3000 \
            --out relatedness

        # Require the expected deliverables before reporting completion.
        for output_file in autosomes.callrate98.pgen autosomes.callrate98.pvar \
            autosomes.callrate98.psam relatedness.kin0; do
            if [[ ! -s "$output_file" ]]; then
                echo "ERROR: an expected output file is missing or empty" >&2
                exit 1
            fi
        done
        echo "Sample call-rate filtering and KING-robust: complete"
    >>>

    output {
        File filtered_pgen = "autosomes.callrate98.pgen"
        File filtered_pvar = "autosomes.callrate98.pvar"
        File filtered_psam = "autosomes.callrate98.psam"
        File king_table = "relatedness.kin0"
        Array[File] plink_logs = ["autosomes.callrate98.log", "relatedness.log"]
    }

    # Part 3. Container and resources
    runtime {
        docker: "ytian820/chip_cad_docker@sha256:fdee1d62eb43ee1e78a2ecf9df55e44a41af6b2b7eb48904ff1db19fa6c57963"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 20 SSD"
    }
}
