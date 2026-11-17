version 1.0

# Purpose: independently merge screened autosomal and X VCFs for the N476 cohort.
# Phenotype: none. Unit: one cohort/sample_set per workflow run.
# Inputs: 01 autosomal and 02 X VCFs, each paired with their tabix indexes.
# Outputs: separate multi-sample bgzip VCFs and tabix indexes for autosomes and X.
# Assumptions: inputs have passed separate checks; GRCh37, chromosomes 1-22 or X,
# sorted records, unique sample names and compatible REF/ALT representations.
# Merge the union of records; absent genotypes stay missing, not reference.
# Do not create new multiallelic records (-m none); no QC or harmonization here.
# Each input array must come from the same Terra sample_set membership order.
# Keep original X genotypes and PAR/non-PAR together; no ploidy recoding.

# Part 1. Cohort inputs and outputs
workflow MergeCohort {
    input {
        Array[File]+ autosomal_vcfs
        Array[File]+ autosomal_vcf_indexes
        Array[File]+ chrX_vcfs
        Array[File]+ chrX_vcf_indexes
    }

    # Independent calls: neither consumes the other's outputs.
    call MergeVcfs as MergeAutosomes {
        input:
            vcfs = autosomal_vcfs,
            vcf_indexes = autosomal_vcf_indexes,
            output_prefix = "cohort.autosomes"
    }

    call MergeVcfs as MergeChrX {
        input:
            vcfs = chrX_vcfs,
            vcf_indexes = chrX_vcf_indexes,
            output_prefix = "cohort.chrX"
    }

    output {
        File merged_autosomal_vcf = MergeAutosomes.merged_vcf
        File merged_autosomal_vcf_index = MergeAutosomes.merged_vcf_index
        File merged_chrX_vcf = MergeChrX.merged_vcf
        File merged_chrX_vcf_index = MergeChrX.merged_vcf_index
    }
}

# Part 2. Stage localized input pairs and merge
task MergeVcfs {
    input {
        Array[File]+ vcfs
        Array[File]+ vcf_indexes
        String output_prefix
    }

    command <<<
        set -euo pipefail

        echo "Purpose: merge screened VCFs for one cohort; output: ~{output_prefix}"
        echo "Input VCF count: ~{length(vcfs)}"
        echo "Merge: union of records; missing genotypes remain missing; -m none"
        echo "Resources: 2 CPUs, 16 GB memory, 100 GB local disk"
        bcftools --version

        # Interpolate File arrays here, after localization. Indexes may have
        # different localized parent directories, so stage each pair together.
        python3 - ~{length(vcfs)} "~{sep='" "' vcfs}" "~{sep='" "' vcf_indexes}" <<'PYTHON'
        import os
        import sys

        n = int(sys.argv[1])
        vcfs = sys.argv[2:2 + n]
        indexes = sys.argv[2 + n:]
        os.mkdir("merge_inputs")
        with open("vcfs.list", "x") as handle:
            for i, (vcf, index) in enumerate(zip(vcfs, indexes, strict=True)):
                staged_vcf = f"merge_inputs/input_{i:04d}.vcf.gz"
                os.symlink(os.path.abspath(vcf), staged_vcf)
                os.symlink(os.path.abspath(index), staged_vcf + ".tbi")
                handle.write(staged_vcf + "\n")
        PYTHON

        bcftools merge --no-version --merge none --threads 1 \
            --file-list vcfs.list \
            --output-type z --output "~{output_prefix}.vcf.gz"
        bcftools index --tbi "~{output_prefix}.vcf.gz"
        echo "Cohort_merge: complete"
    >>>

    output {
        File merged_vcf = "~{output_prefix}.vcf.gz"
        File merged_vcf_index = "~{output_prefix}.vcf.gz.tbi"
    }

    # Part 3. Container and resources
    runtime {
        docker: "ytian820/chip_cad_docker@sha256:fdee1d62eb43ee1e78a2ecf9df55e44a41af6b2b7eb48904ff1db19fa6c57963"
        cpu: 2
        memory: "16 GB"
        disks: "local-disk 100 SSD"
    }
}
