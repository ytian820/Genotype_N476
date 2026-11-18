version 1.0

# Purpose: download PLINK2's 1KG Phase 3 GRCh37 reference for 05b.
# Phenotype: none. Unit: public reference cohort.
# Inputs: none; requires HTTPS access to Dropbox and the EBI archive.
# Outputs: uncompressed PGEN/PVAR/PSAM, official population panel and manifest.
# Sources: https://www.cog-genomics.org/plink/2.0/resources (GRCh37, not hg38).
# Use the full variant set, no INFO annotations, and KING-corrected PSAM.
# No sample/variant filtering or PCA; restrict to autosomes later in 05b.
# SHA256 records downloaded bytes; it is not an upstream checksum comparison.

# Part 1. Workflow interface
workflow Prepare1kgReference {
    call Download1kgReference

    output {
        File reference_pgen = Download1kgReference.reference_pgen
        File reference_pvar = Download1kgReference.reference_pvar
        File reference_psam = Download1kgReference.reference_psam
        File population_panel = Download1kgReference.population_panel
        File download_manifest = Download1kgReference.download_manifest
    }
}

# Part 2. Download and verify the reference bundle
task Download1kgReference {
    command <<<
        set -euo pipefail
        echo "Reference: PLINK2 1KG Phase 3 GRCh37; full variant set; KING-corrected PSAM"
        echo "Scope: download reference trio and official population labels; no filtering"
        echo "Checks: nonempty downloads, HTTP content length when available, Zstandard decompression"
        echo "Resources: 1 CPU, 4 GB memory, 40 GB local disk"
        command -v python3 >/dev/null
        command -v plink2 >/dev/null

        python3 <<'PY'
        import hashlib
        import time
        import urllib.request
        from pathlib import Path

        # URLs linked by the PLINK2 resource page, plus the official EBI panel.
        downloads = [
            ("reference.pgen.zst", "https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1"),
            ("reference.pvar.zst", "https://www.dropbox.com/s/c95n8quqwqww4s0/all_phase3_noannot.pvar.zst?dl=1"),
            ("reference.psam", "https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?dl=1&rlkey=0yyifzj2fb863ddbmsv4jkeq6"),
            ("integrated_call_samples_v3.20130502.ALL.panel", "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"),
        ]
        manifest_rows = ["filename\turl\tbytes\tsha256"]

        def download(name, url):
            partial = Path(name + ".partial")
            for attempt in range(1, 6):
                try:
                    digest = hashlib.sha256()
                    size = 0
                    with urllib.request.urlopen(url, timeout=120) as response:
                        expected_size = response.headers.get("Content-Length")
                        with partial.open("wb") as output:
                            while chunk := response.read(1024 * 1024):
                                output.write(chunk)
                                digest.update(chunk)
                                size += len(chunk)
                    if size == 0 or (expected_size and size != int(expected_size)):
                        raise RuntimeError("Empty or incomplete download")
                    partial.replace(name)
                    manifest_rows.append(f"{name}\t{url}\t{size}\t{digest.hexdigest()}")
                    return
                except Exception:
                    if attempt == 5:
                        raise
                    print(f"Retrying download: {name}; attempt {attempt + 1}/5", flush=True)
                    time.sleep(10)

        # Step 1. Download four files and record their sources and checksums.
        for name, url in downloads:
            download(name, url)
            print(f"Downloaded {name}", flush=True)
        Path("download_manifest.tsv").write_text("\n".join(manifest_rows) + "\n")
        PY

        # Step 2. Decompress using PLINK2; all three files share one prefix.
        plink2 --zst-decompress reference.pgen.zst reference.pgen
        plink2 --zst-decompress reference.pvar.zst reference.pvar
        test -s reference.pgen
        test -s reference.pvar
        echo "Reference preparation complete: PGEN/PVAR/PSAM and population panel"
    >>>

    output {
        File reference_pgen = "reference.pgen"
        File reference_pvar = "reference.pvar"
        File reference_psam = "reference.psam"
        File population_panel = "integrated_call_samples_v3.20130502.ALL.panel"
        File download_manifest = "download_manifest.tsv"
    }

    # Part 3. Container and resources
    runtime {
        docker: "ytian820/chip_cad_docker@sha256:fdee1d62eb43ee1e78a2ecf9df55e44a41af6b2b7eb48904ff1db19fa6c57963"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 40 SSD"
    }
}
