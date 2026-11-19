version 1.0

# Purpose: flag autosomal heterozygosity outliers using PLINK2 F in the 05a cohort.
# Phenotype: none. Unit: genotyped sample, identified by FID/IID.
# Inputs: 05a PGEN/PVAR/PSAM and 05c ancestry assignments.
# Outputs: all-sample table, flagged-sample table and PLINK log.
# Execution only: input checks, plots and summary review are performed manually.
# Parameters: biallelic A/C/G/T SNPs on chromosomes 1-22; hardcalls.
# No extra MAF, missingness or LD filter; PLINK --het skips monomorphic sites.
# Use the F column from --het (method-of-moments inbreeding coefficient).
# Flag F outside mean(F) +/- 3 * sample SD(F): low_F or high_F.
# Low F indicates excess heterozygosity; high F indicates reduced heterozygosity.
# SD uses ddof=1; allele frequencies are estimated from the input cohort.
# Ancestry uses predicted_group for annotation only; no PC adjustment,
# ancestry-specific thresholds or sample removal.

workflow Heterozygosity {
    input {
        File pgen
        File pvar
        File psam
        File ancestry_assignments
    }

    call FlagHeterozygosity {
        input:
            pgen = pgen, pvar = pvar, psam = psam,
            ancestry_assignments = ancestry_assignments
    }

    output {
        File heterozygosity = FlagHeterozygosity.heterozygosity
        File heterozygosity_outliers = FlagHeterozygosity.heterozygosity_outliers
        File plink_log = FlagHeterozygosity.plink_log
    }
}

task FlagHeterozygosity {
    input {
        File pgen
        File pvar
        File psam
        File ancestry_assignments
    }

    command <<<
        set -euo pipefail

        # Part 1. Run configuration (PLINK2 and Python are provided by the image)
        echo "Autosomal biallelic SNP PLINK2 F statistic; all samples retained"
        echo "Flag: F outside whole-cohort mean(F) +/- 3 * sample SD(F); ancestry annotation only"
        echo "No extra MAF/missingness filter, LD pruning or PC adjustment"
        echo "Resources: 1 CPU, 4 GB RAM, 20 GB disk"

        # Part 2. Calculate PLINK2 F (same variant filters as terminal preview)
        plink2 --pgen "~{pgen}" --pvar "~{pvar}" --psam "~{psam}" \
            --autosome --snps-only just-acgt --max-alleles 2 --het \
            --threads 1 --memory 3000 --out cohort

        # Part 3. Flag both tails of F, add ancestry and write tables
        python3 - "~{ancestry_assignments}" <<'PY'
        import csv
        import statistics
        import sys

        def read_table(path):
            with open(path) as handle:
                columns = handle.readline().lstrip('#').split()
                rows = [dict(zip(columns, line.split())) for line in handle if line.strip()]
            for row in rows:
                row.setdefault('FID', '0')
            return rows

        # Keep FID/IID as strings; ancestry only annotates the genotype samples.
        samples = read_table('cohort.het')
        ancestry = {
            (row['FID'], row['IID']): row
            for row in read_table(sys.argv[1])
        }
        f_values = [float(row['F']) for row in samples]
        mean = statistics.mean(f_values)
        sd = statistics.stdev(f_values)
        lower = mean - 3 * sd
        upper = mean + 3 * sd

        for row, value in zip(samples, f_values):
            label = ancestry.get((row['FID'], row['IID']), {})
            row.update(
                f_z_score=(value - mean) / sd if sd > 0 else 0.0,
                low_F=value < lower,
                high_F=value > upper,
                predicted_group=label.get('predicted_group', ''),
                max_probability=label.get('max_probability', '')
            )

        flagged = [row for row in samples if row['low_F'] or row['high_F']]
        for path, rows in [('heterozygosity.tsv', samples),
                           ('heterozygosity_outliers.tsv', flagged)]:
            with open(path, 'w', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(samples[0]), delimiter='\t')
                writer.writeheader()
                writer.writerows(rows)
        PY
    >>>

    output {
        File heterozygosity = "heterozygosity.tsv"
        File heterozygosity_outliers = "heterozygosity_outliers.tsv"
        File plink_log = "cohort.log"
    }

    runtime {
        docker: "ytian820/chip_cad_docker@sha256:fdee1d62eb43ee1e78a2ecf9df55e44a41af6b2b7eb48904ff1db19fa6c57963"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 20 SSD"
    }
}
