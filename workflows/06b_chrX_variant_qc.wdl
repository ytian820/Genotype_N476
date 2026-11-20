version 1.0

# Purpose: chrX variant QC before TOPMed imputation; GRCh37 hardcalls.
# Phenotype: none. Units: genotyped sample (FID/IID) and variant.
# Inputs: complete 03 chrX VCF, 05f final keep list, reviewed 05e sex labels.
# Outputs: filtered trio with sex, pre-filter statistics and PLINK logs.
# Parameters: biallelic A/C/G/T SNPs; geno <= 0.02; MAF >= 0.01;
# HWE 1e-5 midp keep-fewhet (k=0): default joint test on non-PAR X,
# diploid test on PAR. Frequencies and HWE use founders only.
# No parental IDs are imported, so all retained samples are treated as founders.
# Split PAR before clearing male non-PAR heterozygous hardcalls.
# Preserve variant IDs and alleles; no LD pruning or additional sample exclusions.

workflow ChrXVariantQC {
    input {
        File chrX_vcf
        File samples_keep
        File sex_labels
    }

    call FilterChrXVariants {
        input: chrX_vcf = chrX_vcf,
               samples_keep = samples_keep, sex_labels = sex_labels
    }

    output {
        File filtered_pgen = FilterChrXVariants.filtered_pgen
        File filtered_pvar = FilterChrXVariants.filtered_pvar
        File filtered_psam = FilterChrXVariants.filtered_psam
        File variant_missingness = FilterChrXVariants.variant_missingness
        File allele_frequencies = FilterChrXVariants.allele_frequencies
        Array[File] hwe_statistics = FilterChrXVariants.hwe_statistics
        Array[File] plink_logs = FilterChrXVariants.plink_logs
    }
}

task FilterChrXVariants {
    input {
        File chrX_vcf
        File samples_keep
        File sex_labels
    }

    command <<<
        set -euo pipefail

        # Part 1. Prepare sex codes for the final sample set
        Rscript - "~{samples_keep}" "~{sex_labels}" <<'RSCRIPT'
        suppressPackageStartupMessages(library(tidyverse))

        args <- commandArgs(trailingOnly = TRUE)
        keep <- read_tsv(args[1], col_types = cols(.default = col_character())) |>
            rename(FID = `#FID`)
        labels <- read_tsv(args[2], col_types = cols(.default = col_character())) |>
            select(FID, IID, inferred_sex)

        stopifnot(!anyDuplicated(labels[c("FID", "IID")]))
        sex <- left_join(keep, labels, by = c("FID", "IID"))
        if (!all(sex$inferred_sex %in% c("Male", "Female"))) {
            stop("Every retained sample must have a reviewed Male/Female label")
        }
        sex <- sex |>
            transmute(`#FID` = FID, IID, SEX = if_else(inferred_sex == "Male", 1L, 2L))
        write_tsv(sex, "sex.update")
        RSCRIPT

        # Part 2. Import SNPs, update sex and split PAR before clearing invalid calls
        plink2 --vcf "~{chrX_vcf}" --double-id --vcf-half-call missing \
            --keep "~{samples_keep}" --update-sex sex.update --split-par b37 \
            --snps-only just-acgt --min-alleles 2 --max-alleles 2 \
            --make-pgen --threads 1 --memory 3000 --out imported

        plink2 --pfile imported --chr X,PAR1,PAR2 \
            --set-invalid-haploid-missing --make-pgen \
            --threads 1 --memory 3000 --out prepared

        # Part 3. Report statistics after genotype cleanup, before numerical filters
        # PLINK writes separate HWE reports for PAR and non-PAR when present.
        plink2 --pfile prepared --freq --missing variant-only --hardy midp \
            --threads 1 --memory 3000 --out before_filters

        # Part 4. Apply variant filters and write the final trio
        # PLINK uses its internal filter order, regardless of flag order.
        plink2 --pfile prepared \
            --geno 0.02 --maf 0.01 --hwe 1e-5 midp keep-fewhet \
            --make-pgen --threads 1 --memory 3000 --out chrX.qc
    >>>

    output {
        File filtered_pgen = "chrX.qc.pgen"
        File filtered_pvar = "chrX.qc.pvar"
        File filtered_psam = "chrX.qc.psam"
        File variant_missingness = "before_filters.vmiss"
        File allele_frequencies = "before_filters.afreq"
        Array[File] hwe_statistics = glob("before_filters.hardy*")
        Array[File] plink_logs = ["imported.log", "prepared.log", "before_filters.log", "chrX.qc.log"]
    }

    runtime {
        # chip_cad_docker:1.0.5
        docker: "ytian820/chip_cad_docker@sha256:0a3c980cfc1f6bd6f8eb47c3993864976116a09b6f0b974dc7107aabd2619c34"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 20 SSD"
    }
}
