version 1.0

# Purpose: test the CHIP-CAD image in Terra without reading research data.
# Inputs: none. Output: tool versions in the task's stdout log.
workflow TestTools {
    call PrintVersions
}

task PrintVersions {
    command <<<
        set -e
        python3 --version
        bcftools --version
        plink2 --version
    >>>

    runtime {
        docker: "ytian820/chip_cad_docker@sha256:fdee1d62eb43ee1e78a2ecf9df55e44a41af6b2b7eb48904ff1db19fa6c57963"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 10 SSD"
    }
}
