version 1.0

# Purpose: generate the reviewed sample list for downstream variant QC.
# Phenotype: none. Unit: genotyped sample, matched by exact FID/IID.
# Inputs: 05e keep list (already excludes the reviewed heterozygosity outlier)
# and the manual sex-discordance list; first two columns must be FID/IID.
# Both lists may have a FID/IID or #FID/IID header.
# Output: PLINK keep list (#FID/IID), preserving input order.
# Remove only listed sex-discordant samples; no additional QC thresholds.
# Execution only: input linkage and output counts are reviewed manually.

workflow FinalizeSamples {
    input {
        File previous_samples_keep
        File sex_discordant_samples
    }

    call GenerateSampleKeepList {
        input: previous_samples_keep = previous_samples_keep,
               sex_discordant_samples = sex_discordant_samples
    }

    output {
        File final_samples_keep = GenerateSampleKeepList.final_samples_keep
    }
}

task GenerateSampleKeepList {
    input {
        File previous_samples_keep
        File sex_discordant_samples
    }

    command <<<
        set -euo pipefail

        awk '
            BEGIN {
                OFS = "\t"
                print "#FID", "IID"
            }

            # Normalize CRLF line endings before matching FID/IID.
            { sub(/\r$/, "") }

            # Skip headers and blank lines in either file.
            FNR == 1 && ($1 == "#FID" || $1 == "FID") { next }
            NF < 2 { next }

            # First file: remember the FID/IID pairs to exclude.
            FILENAME == ARGV[1] {
                excluded[$1, $2] = 1
                next
            }

            # Second file: keep samples not listed for exclusion.
            {
                fid = $1
                iid = $2
                if (!((fid, iid) in excluded)) {
                    print fid, iid
                }
            }
        ' "~{sex_discordant_samples}" "~{previous_samples_keep}" > final_samples.keep
    >>>

    output {
        File final_samples_keep = "final_samples.keep"
    }

    runtime {
        docker: "ytian820/chip_cad_docker@sha256:fdee1d62eb43ee1e78a2ecf9df55e44a41af6b2b7eb48904ff1db19fa6c57963"
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 10 SSD"
    }
}
