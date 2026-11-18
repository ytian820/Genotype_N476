version 1.0

# Purpose: check 03 cohort merges against the screened single-sample inputs.
# Phenotype: none. Unit: one N476 sample_set; autosomes and X checked separately.
# Inputs: merged and original screened VCFs, each paired with a tabix index.
# Outputs: JSON reports and an overall Boolean; no input VCF is modified.
# Parameters: 476 samples; chromosomes 1-22 or X; sorted, unique biallelic SNPs.
# Compare exact GT strings at original sites (including phase and ploidy).
# At sites absent from a sample's input, its merged GT must be fully missing.
# PASS does not establish biological identity or validate the reference build.
# A completed Terra task is not a QC PASS: inspect all_passed / summary.json.

# Part 1. Independent chromosome checks and combined summary
workflow Check03BasicCheck {
    input {
        Array[File]+ autosomal_vcfs
        Array[File]+ autosomal_vcf_indexes
        Array[File]+ chrX_vcfs
        Array[File]+ chrX_vcf_indexes
        File merged_autosomal_vcf
        File merged_autosomal_vcf_index
        File merged_chrX_vcf
        File merged_chrX_vcf_index
    }

    call CheckMerge as CheckAutosomes {
        input: vcfs = autosomal_vcfs, indexes = autosomal_vcf_indexes,
            merged = merged_autosomal_vcf, merged_index = merged_autosomal_vcf_index,
            group = "autosomes"
    }
    call CheckMerge as CheckChrX {
        input: vcfs = chrX_vcfs, indexes = chrX_vcf_indexes,
            merged = merged_chrX_vcf, merged_index = merged_chrX_vcf_index,
            group = "chrX"
    }
    call SummarizeChecks {
        input: autosomes = CheckAutosomes.report, chrX = CheckChrX.report
    }
    output {
        File autosomal_report = CheckAutosomes.report
        File chrX_report = CheckChrX.report
        File summary = SummarizeChecks.summary
        Boolean all_passed = SummarizeChecks.all_passed
    }
}

# Part 2. Stream inputs, then stream the merged VCF once
# Store site membership as a bit mask, not a cohort-sized genotype matrix.
task CheckMerge {
    input {
        Array[File]+ vcfs
        Array[File]+ indexes
        File merged
        File merged_index
        String group
    }
    command <<<
        set -euo pipefail
        python3 - "~{group}" "~{merged}" "~{merged_index}" ~{length(vcfs)} "~{sep='" "' vcfs}" "~{sep='" "' indexes}" <<'PYTHON'
        import gzip
        import hashlib
        import json
        import os
        import re
        import subprocess
        import sys
        from pathlib import Path

        group, merged, merged_index = sys.argv[1:4]
        n = int(sys.argv[4])
        vcfs, indexes = sys.argv[5:5+n], sys.argv[5+n:]
        expected_samples = 476
        chromosomes = [str(i) for i in range(1, 23)] if group == "autosomes" else ["X"]
        rank = {chrom: i for i, chrom in enumerate(chromosomes)}
        gt_pattern = re.compile(r"[01.](?:[/|][01.])*")
        report = {"group": group, "status": "FAIL", "expected_samples": expected_samples}
        print(f"Purpose: check merged {group}; expected samples: {expected_samples}", flush=True)
        print("Checks: sample names, union of SNPs, exact GT and indexed retrieval", flush=True)

        def require(condition, message):
            if not condition:
                raise ValueError(message)

        def stage(vcf, index, label):
            require(Path(vcf).stat().st_size > 0 and Path(index).stat().st_size > 0,
                    "Empty VCF or index")
            path = f"{label}.vcf.gz"
            os.symlink(os.path.abspath(vcf), path)
            os.symlink(os.path.abspath(index), path + ".tbi")
            return path

        def header(handle):
            for line in handle:
                if line.startswith("##"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                require(fields[:9] == ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"],
                        "Invalid VCF header")
                names = fields[9:]
                require(names and all(names) and len(set(names)) == len(names),
                        "Empty or duplicate sample names")
                return names
            raise ValueError("Missing VCF header")

        def records(handle, sample_count):
            previous = (-1, 0)
            for line in handle:
                fields = line.rstrip("\r\n").split("\t")
                require(len(fields) == 9 + sample_count, "Invalid VCF column count")
                chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
                require(chrom in rank and pos.isdecimal() and int(pos) > 0, "Invalid chromosome or position")
                coordinate = (rank[chrom], int(pos))
                require(coordinate > previous, "Unsorted or duplicate coordinate")
                previous = coordinate
                require(ref in ("A", "C", "G", "T") and alt in ("A", "C", "G", "T") and ref != alt,
                        "Not a distinct biallelic SNP")
                require(fields[8] == "GT", "Expected GT-only FORMAT")
                gts = fields[9:]
                require(all(gt_pattern.fullmatch(gt) for gt in gts), "Invalid GT")
                yield (chrom, str(int(pos)), ref, alt), gts

        def site_bytes(key):
            return ("\t".join(key) + "\n").encode()

        def update_gt(digest, key_bytes, gt):
            digest.update(key_bytes)
            digest.update((gt + "\n").encode())

        def check_index(path, sequential_digest):
            # Query all expected chromosomes through the supplied index.
            digest = hashlib.sha256()
            with subprocess.Popen(
                ["bcftools", "query", "-r", ",".join(chromosomes),
                 "-f", "%CHROM\t%POS\t%REF\t%ALT\n", path],
                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
            ) as process:
                for line in process.stdout:
                    digest.update(line)
                require(process.wait() == 0, "Indexed query failed")
            require(digest.digest() == sequential_digest.digest(), "Index retrieval differs from full scan")

        def run():
            require(n == expected_samples and len(indexes) == n, "Expected 476 VCF/index pairs")
            require(len(set(vcfs)) == n and len(set(indexes)) == n, "Repeated input VCF or index path")
            names, expected_hashes = [], []
            membership = {}
            # Each bit identifies the input sample that contains this exact site.
            for i, (vcf, index) in enumerate(zip(vcfs, indexes)):
                path = stage(vcf, index, f"input_{i}")
                gt_hash, site_hash = hashlib.sha256(), hashlib.sha256()
                count = 0
                with gzip.open(path, "rt") as handle:
                    sample_names = header(handle)
                    require(len(sample_names) == 1, "Input is not single-sample")
                    require(sample_names[0] not in names, "Duplicate input sample name")
                    names.append(sample_names[0])
                    for key, gts in records(handle, 1):
                        membership[key] = membership.get(key, 0) | (1 << i)
                        key_bytes = site_bytes(key)
                        site_hash.update(key_bytes)
                        update_gt(gt_hash, key_bytes, gts[0])
                        count += 1
                require(count > 0, "Empty input variant set")
                check_index(path, site_hash)
                expected_hashes.append(gt_hash.digest())
                if (i + 1) % 25 == 0 or i + 1 == n:
                    print(f"Checked input files: {i + 1}/{n}", flush=True)

            report["input_union_sites"] = len(membership)
            observed_hashes = [hashlib.sha256() for _ in names]
            path = stage(merged, merged_index, "merged")
            site_hash = hashlib.sha256()
            count = 0
            with gzip.open(path, "rt") as handle:
                merged_names = header(handle)
                require(len(merged_names) == n and set(merged_names) == set(names), "Merged sample names differ from inputs")
                # Match by internal VCF name; Terra row IDs need not use this spelling.
                input_order = {name: i for i, name in enumerate(names)}
                columns = [input_order[name] for name in merged_names]
                for key, gts in records(handle, n):
                    require(key in membership, "Merged site not in input union (or alleles differ)")
                    mask = membership.pop(key)
                    key_bytes = site_bytes(key)
                    site_hash.update(key_bytes)
                    for column, gt in zip(columns, gts):
                        if mask & (1 << column):
                            update_gt(observed_hashes[column], key_bytes, gt)
                        else:
                            require(all(allele == "." for allele in re.split(r"[/|]", gt)),
                                    "Nonmissing GT at a site absent from that sample's input")
                    count += 1
            require(not membership, "Input sites missing from merged VCF")
            check_index(path, site_hash)
            mismatches = sum(observed.digest() != expected for observed, expected in zip(observed_hashes, expected_hashes))
            report.update(merged_samples=len(merged_names), merged_sites=count,
                          genotype_mismatched_samples=mismatches, sample_names=merged_names)
            require(mismatches == 0, "GT changed during merge")
            report["status"] = "PASS"

        try:
            run()
        except ValueError as error:
            report["error"] = str(error)
        except Exception as error:
            # Avoid printing sensitive paths or sample/genotype values.
            report["error"] = type(error).__name__ + ": check files, tools and runtime"
        with open("report.json", "x") as handle:
            json.dump(report, handle, indent=2)
        print(f"{group}: {report['status']}", flush=True)
        PYTHON
    >>>
    output { File report = "report.json" }
    runtime {
        docker: "ytian820/chip_cad_docker@sha256:fdee1d62eb43ee1e78a2ecf9df55e44a41af6b2b7eb48904ff1db19fa6c57963"
        cpu: 2
        memory: "16 GB"
        disks: "local-disk 100 SSD"
    }
}

# Part 3. Require both reports and matching autosomal/X sample order
task SummarizeChecks {
    input { File autosomes File chrX }
    command <<<
        set -euo pipefail
        python3 - "~{autosomes}" "~{chrX}" <<'PYTHON'
        import json
        import sys
        with open(sys.argv[1]) as handle:
            autosomes = json.load(handle)
        with open(sys.argv[2]) as handle:
            chrx = json.load(handle)
        same_order = bool(autosomes.get("sample_names")) and autosomes.get("sample_names") == chrx.get("sample_names")
        passed = autosomes["status"] == chrx["status"] == "PASS" and same_order
        result = {"overall": "PASS" if passed else "FAIL",
                  "autosomes": autosomes["status"], "chrX": chrx["status"],
                  "same_sample_order": same_order,
                  "scope": "File/index, sample names, SNP union and exact GT; no biological identity or reference-build validation"}
        with open("summary.json", "x") as handle:
            json.dump(result, handle, indent=2)
        with open("passed.txt", "x") as handle:
            handle.write("true\n" if passed else "false\n")
        print(json.dumps(result))
        PYTHON
    >>>
    output {
        File summary = "summary.json"
        Boolean all_passed = read_boolean("passed.txt")
    }
    runtime {
        docker: "ytian820/chip_cad_docker@sha256:fdee1d62eb43ee1e78a2ecf9df55e44a41af6b2b7eb48904ff1db19fa6c57963"
        cpu: 1
        memory: "1 GB"
        disks: "local-disk 10 SSD"
    }
}
