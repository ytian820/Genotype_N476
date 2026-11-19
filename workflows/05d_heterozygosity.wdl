version 1.0

# Purpose: flag autosomal heterozygosity outliers using PLINK2 F in the 05a cohort.
# Phenotype: none. Unit: genotyped sample, identified by FID/IID.
# Inputs: 05a PGEN/PVAR/PSAM and 05c ancestry assignments.
# Outputs: all-sample table, flagged-sample table, passing-sample list, retained SNP list and PLINK logs.
# Execution only: input checks, plots and summary review are performed manually.
# Parameters: biallelic A/C/G/T SNPs on chromosomes 1-22; hardcalls.
# MAF >= 0.05, GRCh37 long-range LD exclusion (same regions as 05b),
# then LD pruning (500 kb windows, step 1, r2 <= 0.2).
# No extra missingness filter; calculate F only on the retained SNPs.
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
        File samples_keep = FlagHeterozygosity.samples_keep
        File pruned_variants = FlagHeterozygosity.pruned_variants
        File pruning_log = FlagHeterozygosity.pruning_log
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
        echo "Variants: MAF >= 0.05; GRCh37 long-range LD exclusion; LD pruning 500kb 1 0.2"
        echo "No extra missingness filter, PC adjustment or ancestry-specific thresholds"
        echo "Resources: 1 CPU, 4 GB RAM, 20 GB disk"

        # Part 2. Select common, LD-pruned SNPs and calculate PLINK2 F
        ## Step 1. Apply MAF and long-range LD filters before LD pruning.
        # Same frozen GRCh37 ranges as 05b, 1-based inclusive:
        # https://github.com/meyer-lab-cshl/plinkQC/blob/master/inst/extdata/high-LD-regions-hg19-GRCh37.txt
        cat > long_ld.txt <<'REGIONS'
        1 48227412 52227412 1
        2 86146488 101133568 2
        2 89563673 89564169 3
        2 98080299 98080324 4
        2 183291754 190291755 5
        3 47524995 50024996 6
        3 83417309 86917310 7
        5 44464242 50464243 8
        5 128972100 131972101 9
        6 25392020 33392022 10
        6 56892040 63942041 11
        6 139958306 142458307 12
        7 55032505 66362565 13
        8 7962589 11962591 14
        8 42880842 49837447 15
        8 111930823 114930824 16
        10 36959993 43679994 17
        11 87860351 90860352 18
        12 33108732 41713733 19
        20 32536338 35066586 20
        REGIONS

        plink2 --pgen "~{pgen}" --pvar "~{pvar}" --psam "~{psam}" \
            --autosome --snps-only just-acgt --max-alleles 2 \
            --maf 0.05 --exclude range long_ld.txt \
            --indep-pairwise 500kb 1 0.2 \
            --threads 1 --memory 3000 --out pruning

        ## Step 2. Use only the retained SNPs; pruning does not filter this automatically.
        plink2 --pgen "~{pgen}" --pvar "~{pvar}" --psam "~{psam}" \
            --extract pruning.prune.in --het \
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

        # PLINK --keep list: samples within mean(F) +/- 3 SD, including boundaries.
        with open('samples.keep', 'w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['#FID', 'IID'])
            for row, value in zip(samples, f_values):
                if lower <= value <= upper:
                    writer.writerow([row['FID'], row['IID']])
        PY
    >>>

    output {
        File heterozygosity = "heterozygosity.tsv"
        File heterozygosity_outliers = "heterozygosity_outliers.tsv"
        File samples_keep = "samples.keep"
        File pruned_variants = "pruning.prune.in"
        File pruning_log = "pruning.log"
        File plink_log = "cohort.log"
    }

    runtime {
        docker: "ytian820/chip_cad_docker@sha256:fdee1d62eb43ee1e78a2ecf9df55e44a41af6b2b7eb48904ff1db19fa6c57963"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 20 SSD"
    }
}
