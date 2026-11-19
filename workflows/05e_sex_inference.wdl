version 1.0

# Purpose: infer sex from non-PAR chrX for manual review.
# Phenotype: none. Unit: genotyped sample; preserve VCF IDs as FID and IID.
# Inputs: 03 merged GRCh37 chrX VCF and tabix index; 05a PSAM (476 samples).
# Outputs: 475-sample keep list, X_F and sex labels, retained SNPs and logs.
# Exclude only the manually reviewed 05d extreme-low-F sample; retain high-F EAS.
# Parameters: biallelic A/C/G/T SNPs; geno 0.01; MAF 0.05;
# LD pruning 500kb 1 0.2; no HWE filter or ancestry-specific frequencies.
# Explicit PLINK sex thresholds: max-female-xf=0.2, min-male-xf=0.8.
# No reported-sex input, genotype erasure, sex updates or further sample removal.
# Labels are preliminary: Male / Female / Uncertain; review X_F manually.

workflow SexInference {
    input {
        File chrX_vcf
        File chrX_vcf_index
        File psam
    }

    call InferSex {
        input: chrX_vcf = chrX_vcf, chrX_vcf_index = chrX_vcf_index,
               psam = psam
    }

    output {
        File samples_keep = InferSex.samples_keep
        File sex_labels = InferSex.sex_labels
        File pruned_variants = InferSex.pruned_variants
        Array[File] plink_logs = InferSex.plink_logs
    }
}

task InferSex {
    input {
        File chrX_vcf
        File chrX_vcf_index
        File psam
    }

    command <<<
        set -euo pipefail

        # Part 1. Generate the keep list from the 05a PSAM
        echo "Sex inference: GRCh37 non-PAR X; 476 samples from 05a PSAM, one manual exclusion"
        echo "Variants: biallelic A/C/G/T; geno 0.01; MAF 0.05; no HWE filter"
        echo "LD pruning: 500kb 1 0.2; sex thresholds: female 0.2, male 0.8"
        echo "Outputs: 475-sample keep list and X_F/sex labels; manual review only"
        echo "Resources: 1 CPU, 4 GB RAM, 20 GB disk"

        ln -s "~{chrX_vcf}" input.vcf.gz
        ln -s "~{chrX_vcf_index}" input.vcf.gz.tbi

        python3 - "~{psam}" <<'PY'
        import csv
        import sys

        excluded = '209596310147_R08C01'
        with open(sys.argv[1]) as handle:
            columns = handle.readline().lstrip('#').split()
            rows = [dict(zip(columns, line.split())) for line in handle if line.strip()]
        ids = [(row['FID'], row['IID']) for row in rows]
        if len(ids) != 476 or len(set(ids)) != 476:
            raise ValueError('Expected exactly 476 unique FID/IID pairs in the 05a PSAM')
        retained = [(fid, iid) for fid, iid in ids if iid != excluded]
        if len(retained) != 475:
            raise ValueError('Expected exactly one excluded sample and 475 retained samples')

        with open('samples.keep', 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['#FID', 'IID'])
            writer.writerows(retained)
        PY

        # Part 2. Import chrX, split PAR, then filter and prune non-PAR SNPs
        # Sex is intentionally unknown during inference. Preserve GT; do not erase
        # heterozygous calls based on a presumed sex. --lax-chrx-import permits this.
        # Unique allele-based IDs allow LD pruning; discard duplicate records.
        plink2 --vcf input.vcf.gz --double-id --vcf-half-call missing \
            --lax-chrx-import --split-par b37 --keep samples.keep \
            --snps-only just-acgt --min-alleles 2 --max-alleles 2 \
            --set-all-var-ids '@:#:$r:$a' --rm-dup exclude-all \
            --make-pgen --threads 1 --memory 3000 --out chrX

        # PAR1/PAR2 now have separate chromosome codes and are excluded by --chr X.
        plink2 --pfile chrX --chr X --geno 0.01 --maf 0.05 \
            --indep-pairwise 500kb 1 0.2 \
            --threads 1 --memory 3000 --out pruning

        # Pruning only writes a list; explicitly extract it for sex inference.
        plink2 --pfile chrX --chr X --extract pruning.prune.in \
            --check-sex max-female-xf=0.2 min-male-xf=0.8 cols=fid,xf \
            --threads 1 --memory 3000 --out inferred

        # Part 3. Export labels without comparing to reported sex
        python3 - <<'PY'
        import csv
        import math

        with open('samples.keep') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            expected = {(row['#FID'], row['IID']) for row in reader}
        with open('inferred.sexcheck') as handle:
            columns = handle.readline().lstrip('#').split()
            rows = [dict(zip(columns, line.split())) for line in handle if line.strip()]
        observed = {(row['FID'], row['IID']) for row in rows}
        if len(rows) != 475 or observed != expected:
            raise ValueError('Sex inference must contain exactly the 475 retained samples')

        labels = {'1': 'Male', '2': 'Female', 'NA': 'Uncertain', '0': 'Uncertain'}
        with open('sex_labels.tsv', 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['FID', 'IID', 'X_F', 'inferred_sex'])
            for row in rows:
                value = float(row['F']) if row['F'] != 'NA' else math.nan
                label = labels[row['SNPSEX']] if math.isfinite(value) else 'Uncertain'
                writer.writerow([row['FID'], row['IID'], row['F'], label])
        PY
    >>>

    output {
        File samples_keep = "samples.keep"
        File sex_labels = "sex_labels.tsv"
        File pruned_variants = "pruning.prune.in"
        Array[File] plink_logs = ["chrX.log", "pruning.log", "inferred.log"]
    }

    runtime {
        docker: "ytian820/chip_cad_docker@sha256:fdee1d62eb43ee1e78a2ecf9df55e44a41af6b2b7eb48904ff1db19fa6c57963"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 20 SSD"
    }
}
