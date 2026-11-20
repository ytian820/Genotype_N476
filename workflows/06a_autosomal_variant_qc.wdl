version 1.0

# Purpose: autosomal variant QC before TOPMed imputation; GRCh37 hardcalls.
# Phenotype: none. Units: genotyped sample (FID/IID) and variant.
# Inputs: complete 05a PGEN/PVAR/PSAM; samples_keep maps to this.05f_final_samples_keep.
# Outputs: filtered trio, pre-filter variant statistics and PLINK logs.
# Execution only: sample matching, duplicate sites and output review are checked manually.
# Parameters: biallelic A/C/G/T SNPs; geno <= 0.01; MAF >= 0.01;
# HWE 1e-5 keep-fewhet (default k=0). Frequencies use the retained cohort.
# No additional sample exclusions, LD pruning, long-range LD exclusion or strand flips.
# Preserve variant IDs and alleles. Reference alignment and VCF export belong to 06b.
# HWE is pooled across ancestries; keep-fewhet retains heterozygote-deficient sites.

workflow AutosomalVariantQC {
    input {
        File pgen
        File pvar
        File psam
        File samples_keep
    }

    call FilterAutosomalVariants {
        input: pgen = pgen, pvar = pvar, psam = psam, samples_keep = samples_keep
    }

    output {
        File filtered_pgen = FilterAutosomalVariants.filtered_pgen
        File filtered_pvar = FilterAutosomalVariants.filtered_pvar
        File filtered_psam = FilterAutosomalVariants.filtered_psam
        File variant_missingness = FilterAutosomalVariants.variant_missingness
        File allele_frequencies = FilterAutosomalVariants.allele_frequencies
        File hwe_statistics = FilterAutosomalVariants.hwe_statistics
        Array[File] plink_logs = FilterAutosomalVariants.plink_logs
    }
}

task FilterAutosomalVariants {
    input {
        File pgen
        File pvar
        File psam
        File samples_keep
    }

    command <<<
        set -euo pipefail

        # Part 1. Run configuration
        echo "Autosomal variant QC: GRCh37; apply the final reviewed 05f keep list"
        echo "Filters: biallelic A/C/G/T SNPs; geno 0.01; MAF 0.01; HWE 1e-5 keep-fewhet (k=0)"
        echo "Statistics: retained cohort; no additional sample filtering or LD pruning"
        echo "Resources: 1 CPU, 4 GB RAM, 20 GB disk"

        # Part 2. Select samples and SNPs; report statistics before numerical filters
        plink2 --pgen "~{pgen}" --pvar "~{pvar}" --psam "~{psam}" \
            --keep "~{samples_keep}" --autosome --snps-only just-acgt \
            --min-alleles 2 --max-alleles 2 --make-pgen \
            --threads 1 --memory 3000 --out selected

        plink2 --pfile selected --nonfounders --freq --missing variant-only --hardy \
            --threads 1 --memory 3000 --out before_filters

        # Part 3. Apply all variant filters and write the final trio
        # PLINK uses its internal filter order, regardless of flag order.
        plink2 --pfile selected --nonfounders \
            --geno 0.01 --maf 0.01 --hwe 1e-5 keep-fewhet \
            --make-pgen --threads 1 --memory 3000 --out autosomes.qc
    >>>

    output {
        File filtered_pgen = "autosomes.qc.pgen"
        File filtered_pvar = "autosomes.qc.pvar"
        File filtered_psam = "autosomes.qc.psam"
        File variant_missingness = "before_filters.vmiss"
        File allele_frequencies = "before_filters.afreq"
        File hwe_statistics = "before_filters.hardy"
        Array[File] plink_logs = ["selected.log", "before_filters.log", "autosomes.qc.log"]
    }

    runtime {
        docker: "ytian820/chip_cad_docker@sha256:fdee1d62eb43ee1e78a2ecf9df55e44a41af6b2b7eb48904ff1db19fa6c57963"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 20 SSD"
    }
}
