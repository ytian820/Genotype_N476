version 1.0

# Purpose: screen chromosome X in one N476 VCF before cohort-level QC.
# Phenotype: none. Unit: one single-sample GRCh37 VCF per workflow run.
# Input: gzip-compressed VCF; Terra localizes it (no gcloud or input index needed).
# Outputs: GT-only chromosome X VCF and tabix index; screening counts in stdout.
# Rules: chromosome named X (not chrX); exclude FILTER tags
# DUPE/TRIALLELIC/ZEROED_OUT_ASSAY; distinct single-base A/C/G/T REF and ALT;
# exclude ALL records at duplicated candidate coordinates within this sample.
# Preserve original sample ID, record order, allele order and GT; remove INFO.
# Keep PAR and non-PAR together; no sex inference or ploidy recoding.
# Run once per selected Terra sample row, mapping vcf to the row's VCF file.
# Fixed output names are isolated in each task's execution directory.

# Part 1. Workflow inputs and outputs
workflow ScreenChrX {
    input {
        File vcf
    }

    call ScreenChrXVcf {
        input: vcf = vcf
    }

    output {
        File screened_vcf = ScreenChrXVcf.screened_vcf
        File screened_vcf_index = ScreenChrXVcf.screened_vcf_index
    }
}

# Part 2. Screen one localized VCF
task ScreenChrXVcf {
    input {
        File vcf
    }

    command <<<
        set -euo pipefail

        echo "Purpose: chromosome X screening of one single-sample GRCh37 VCF"
        echo "Chromosome: X; distinct single-base A/C/G/T REF and ALT"
        echo "Excluded FILTER tags: DUPE,TRIALLELIC,ZEROED_OUT_ASSAY"
        echo "Duplicate candidate coordinates: exclude all records within sample"
        echo "Sample fields: original GT only; INFO removed"
        echo "PAR and non-PAR: both retained; no ploidy recoding"
        echo "Resources: 1 CPU, 4 GB memory, 20 GB local disk"
        python3 --version
        bcftools --version

        # Cromwell substitutes the localized File path in the command.
        python3 - "~{vcf}" <<'PYTHON'
        import gzip
        import json
        import sys
        from collections import Counter

        input_path = sys.argv[1]
        excluded_filters = {"DUPE", "TRIALLELIC", "ZEROED_OUT_ASSAY"}
        bases = {"A", "C", "G", "T"}


        def screen_chrX(input_path, output_path):
            counts = Counter(input=0, excluded_non_X=0,
                             excluded_FILTER=0, excluded_non_SNP=0,
                             excluded_duplicate_coordinate_records=0, retained=0)
            coordinate_counts = Counter()
            candidates = []
            headers = []

            # Step 1. Keep the needed headers and collect candidate SNPs.
            with gzip.open(input_path, "rt") as handle:
                for line in handle:
                    if line.startswith("##INFO="):
                        continue
                    if line.startswith("##FORMAT=") and not line.startswith("##FORMAT=<ID=GT,"):
                        continue
                    if line.startswith("#CHROM") and len(line.rstrip("\r\n").split("\t")) != 10:
                        raise ValueError("Expected a single-sample VCF with exactly 10 columns")
                    if line.startswith("#"):
                        headers.append(line)
                        continue

                    fields = line.rstrip("\r\n").split("\t")
                    counts["input"] += 1
                    chrom = fields[0]
                    if chrom != "X":
                        counts["excluded_non_X"] += 1
                        continue
                    if excluded_filters.intersection(fields[6].split(";")):
                        counts["excluded_FILTER"] += 1
                        continue

                    ref, alt = fields[3].upper(), fields[4].upper()
                    if ref not in bases or alt not in bases or ref == alt:
                        counts["excluded_non_SNP"] += 1
                        continue

                    # Preserve original GT; replace INFO with '.' and FORMAT with 'GT'.
                    coordinate = (chrom, int(fields[1]))
                    gt_index = fields[8].split(":").index("GT")
                    gt = fields[9].split(":")[gt_index]
                    record = "\t".join(fields[:7] + [".", "GT", gt]) + "\n"
                    candidates.append((coordinate, record))
                    coordinate_counts[coordinate] += 1

            # Step 2. Write only records whose candidate coordinate occurs once.
            with open(output_path, "x") as handle:
                handle.writelines(headers)
                for coordinate, record in candidates:
                    if coordinate_counts[coordinate] == 1:
                        handle.write(record)
                        counts["retained"] += 1
                    else:
                        counts["excluded_duplicate_coordinate_records"] += 1
            return counts


        counts = screen_chrX(input_path, "screened.vcf")
        print("Screening_counts: " + json.dumps(counts), flush=True)
        PYTHON

        # Part 3. Compress, index and verify output files
        bcftools view --no-version -Oz -o chrX.screened.vcf.gz screened.vcf
        bcftools index --tbi chrX.screened.vcf.gz
        test -s chrX.screened.vcf.gz
        test -s chrX.screened.vcf.gz.tbi
        echo "ChrX_screening: complete"
    >>>

    output {
        File screened_vcf = "chrX.screened.vcf.gz"
        File screened_vcf_index = "chrX.screened.vcf.gz.tbi"
    }

    # Part 4. Container and resources
    runtime {
        docker: "ytian820/chip_cad_docker@sha256:fdee1d62eb43ee1e78a2ecf9df55e44a41af6b2b7eb48904ff1db19fa6c57963"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 20 SSD"
    }
}
