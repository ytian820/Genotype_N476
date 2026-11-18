version 1.0

# Purpose: remove autosomal variants with genotype missingness > 5%.
# Phenotype: none. Unit: one cohort/sample_set; all input samples are retained.
# Input: checked, sorted, GT-only autosomal VCF from 03 (chromosomes 1-22).
# Outputs: PGEN/PVAR/PSAM, filtered bgzip VCF and tabix index.
# Half-missing GTs are imported as missing (--vcf-half-call missing).
# Sample IDs are imported intact as IID and exported with id-paste=iid.
# PLINK2 re-exports the VCF; original text/annotations are not guaranteed unchanged.
# No MAF, HWE, LD or additional chromosome filter is applied.

# Part 1. Workflow interface
workflow FilterAutosomes {
    input {
        File vcf
    }

    call FilterMissingness {
        input: vcf = vcf
    }

    output {
        File filtered_pgen = FilterMissingness.filtered_pgen
        File filtered_pvar = FilterMissingness.filtered_pvar
        File filtered_psam = FilterMissingness.filtered_psam
        File filtered_vcf = FilterMissingness.filtered_vcf
        File filtered_vcf_index = FilterMissingness.filtered_vcf_index
    }
}

# Part 2. Filter with PLINK2, then export the retained dataset
task FilterMissingness {
    input {
        File vcf
    }

    command <<<
        set -euo pipefail
        echo "Purpose: autosomal variant filtering; retain all samples"
        echo "Threshold: missing GT calls / all samples <= 0.05; half-calls count as missing"
        echo "Resources: 1 CPU, 4 GB memory, 20 GB local disk"

        # Import and filter variants only; no --mind or sample exclusions.
        plink2 --vcf "~{vcf}" \
            --double-id --vcf-half-call missing \
            --geno 0.05 \
            --make-pgen \
            --threads 1 --memory 3000 \
            --out autosomes.missingness05

        # Export from the filtered PGEN files; retain the original IID as VCF ID.
        plink2 --pfile autosomes.missingness05 \
            --export vcf bgz id-paste=iid \
            --threads 1 --memory 3000 \
            --out autosomes.missingness05.export

        bcftools index --tbi autosomes.missingness05.export.vcf.gz
    >>>

    output {
        File filtered_pgen = "autosomes.missingness05.pgen"
        File filtered_pvar = "autosomes.missingness05.pvar"
        File filtered_psam = "autosomes.missingness05.psam"
        File filtered_vcf = "autosomes.missingness05.export.vcf.gz"
        File filtered_vcf_index = "autosomes.missingness05.export.vcf.gz.tbi"
    }

    # Part 3. Container and resources
    runtime {
        docker: "ytian820/chip_cad_docker@sha256:fdee1d62eb43ee1e78a2ecf9df55e44a41af6b2b7eb48904ff1db19fa6c57963"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 20 SSD"
    }
}
