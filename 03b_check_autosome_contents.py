#!/usr/bin/env python3
"""User-run, read-only checks of 476 screened autosomal VCFs before merging.

Phenotype: none. Unit: one sample table row and its single-sample VCF.
Inputs: a fresh Terra sample TSV export, 03 workflow_checks.tsv, cloud VCF/index.
Outputs: content_checks.tsv and summary.txt in a new report directory.
Rules: valid single-sample GT-only autosomal SNPs, sorted unique coordinates,
actual count matching stdout, unique sample names, no REF/ALT conflicts.
Index access must reproduce the full ordered CHROM/POS/REF/ALT digest.
Different site sets are reported, not automatically rejected. ID differences
require review; exact string agreement is not proof of biological identity.
No reference-build validation, genotype concordance or cohort QC is performed.
Run manually in Terra Terminal with Python 3, gcloud and bcftools available.
Downloads one VCF/index pair at a time into a temporary directory, then removes
only those temporary copies. Cloud source files and existing reports are untouched.
"""

import csv
import gzip
import hashlib
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


# Part 1. Fixed configuration
PROJECT_DIR = Path("/home/jupyter/CHIP_CAD_N476")
# Export DATA -> sample as TSV after 01 outputs have been written to the table.
SAMPLE_TABLE = PROJECT_DIR / "data/sample_for_autosome_check.tsv"
COUNT_REPORT = PROJECT_DIR / "result/03_check_autosomes/workflow_checks.tsv"
REPORT_DIR = PROJECT_DIR / "result/03b_check_autosome_contents"
VCF_COLUMN = "01_screen_autosomes_vcf"
INDEX_COLUMN = "01_screen_autosomes_vcf_index"
EXPECTED_SAMPLES = 476
AUTOSOMES = {str(i) for i in range(1, 23)}
REGIONS = ",".join(str(i) for i in range(1, 23))
EXCLUDED_FILTERS = {"DUPE", "TRIALLELIC", "ZEROED_OUT_ASSAY"}


# Part 2. Helpers: manifest, full scan, index access
def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def prepare_manifest():
    columns, samples = read_tsv(SAMPLE_TABLE)
    id_columns = [key for key in ("entity:sample_id", "sample_id") if key in columns]
    if len(id_columns) != 1 or any(key not in columns for key in (VCF_COLUMN, INDEX_COLUMN)):
        raise RuntimeError("Sample TSV needs one sample ID column and both 01 output columns.")
    if len(samples) != EXPECTED_SAMPLES:
        raise RuntimeError("Sample TSV must contain exactly 476 rows.")
    id_column = id_columns[0]
    for key in (id_column, VCF_COLUMN, INDEX_COLUMN):
        values = [row[key] for row in samples]
        if any(not value or value.strip() != value for value in values) or len(set(values)) != len(values):
            raise RuntimeError("Sample IDs and output paths must be populated and unique.")
    for row in samples:
        for key in (VCF_COLUMN, INDEX_COLUMN):
            if not row[key].startswith("gs://") or any(c in row[key] for c in "*?[]#\n\r"):
                raise RuntimeError("Expected exact, unversioned GCS output paths.")

    columns, checks = read_tsv(COUNT_REPORT)
    if any(key not in columns for key in ("workflow_id", "status", "retained", "vcf_uri", "index_uri")):
        raise RuntimeError("Unexpected 03 count-report schema.")
    by_vcf = {row["vcf_uri"]: row for row in checks}
    if len(checks) != EXPECTED_SAMPLES or len(by_vcf) != EXPECTED_SAMPLES:
        raise RuntimeError("03 report must contain 476 unique VCF paths.")
    if set(by_vcf) != {row[VCF_COLUMN] for row in samples}:
        raise RuntimeError("Sample table paths do not match this 03 submission report.")
    for row in samples:
        check = by_vcf[row[VCF_COLUMN]]
        if check["status"] != "PASS" or check["index_uri"] != row[INDEX_COLUMN]:
            raise RuntimeError("03 checks must pass and each index path must match the table.")
        if int(check["retained"]) <= 0:
            raise RuntimeError("03 report contains invalid retained counts.")
    return id_column, samples, by_vcf


def scan_lines(lines):
    """Return counts, identity, ordered site digest and coordinate/allele map."""
    digest = hashlib.sha256()
    sites = {}
    issues = set()
    names = []
    header_count = 0
    records = 0
    previous = None
    for line in lines:
        if line.startswith("#CHROM\t"):
            header_count += 1
            names = line.rstrip("\r\n").split("\t")[9:]
        if line.startswith("#"):
            continue
        records += 1
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) != 10:
            issues.add("invalid_record_columns")
            continue
        chrom, pos, _, ref, alt = fields[:5]
        if chrom not in AUTOSOMES or not pos.isdecimal() or int(pos) <= 0:
            issues.add("invalid_autosomal_coordinate")
            continue
        coordinate = (int(chrom), int(pos))
        if previous is not None and coordinate < previous:
            issues.add("unsorted_coordinates")
        previous = coordinate
        if coordinate in sites:
            issues.add("duplicate_coordinates")
        sites[coordinate] = (ref, alt)
        if ref not in {"A", "C", "G", "T"} or alt not in {"A", "C", "G", "T"} or ref == alt:
            issues.add("invalid_biallelic_snp")
        if fields[7] != "." or fields[8] != "GT":
            issues.add("unexpected_info_or_format")
        if not re.fullmatch(r"[01.](?:[/|][01.])*", fields[9]):
            issues.add("invalid_gt_encoding")
        if EXCLUDED_FILTERS.intersection(fields[6].split(";")):
            issues.add("excluded_filter_still_present")
        digest.update(f"{chrom}\t{pos}\t{ref}\t{alt}\n".encode())
    if header_count != 1 or len(names) != 1 or not names[0]:
        issues.add("expected_one_sample_header")
    if records == 0:
        issues.add("zero_records")
    return names, records, sites, digest.hexdigest(), issues


def download(uri, destination):
    result = subprocess.run(
        ["gcloud", "storage", "cp", uri, str(destination)],
        capture_output=True, text=True, timeout=600,
    )
    if result.returncode:
        raise RuntimeError("Cloud download failed; check account, permissions and network.")


# Part 3. Check one sample at a time
def main():
    print("Purpose: autosomal content checks; expected samples: 476", flush=True)
    print(f"Sample table: {SAMPLE_TABLE}")
    print(f"Count report: {COUNT_REPORT}")
    print(f"Report directory: {REPORT_DIR}")
    print("Checks: IDs, records, order, alleles and indexed retrieval; no genotype comparison.")
    for tool in ("gcloud", "bcftools"):
        if not shutil.which(tool):
            raise RuntimeError(f"Required executable unavailable: {tool}")
    if REPORT_DIR.exists():
        raise RuntimeError("Report directory already exists; nothing overwritten.")
    id_column, samples, by_vcf = prepare_manifest()
    rows = []
    reference_sites = None
    # Track all coordinates, including those absent from the first sample.
    cohort_alleles = {}
    sample_name_rows = {}
    with tempfile.TemporaryDirectory(prefix="autosome_content_check_") as tmp:
        vcf = Path(tmp) / "screened.vcf.gz"
        index = Path(str(vcf) + ".tbi")
        for number, sample in enumerate(samples, start=1):
            check = by_vcf[sample[VCF_COLUMN]]
            download(sample[VCF_COLUMN], vcf)
            download(sample[INDEX_COLUMN], index)
            with gzip.open(vcf, "rt") as handle:
                names, records, sites, digest, issues = scan_lines(handle)
            expected = int(check["retained"])
            if records != expected:
                issues.add("record_count_differs_from_stdout")
            name = names[0] if len(names) == 1 else ""
            review = "sample_id_string_differs" if name and name != sample[id_column] else ""
            if name:
                sample_name_rows.setdefault(name, []).append(len(rows))

            # Query all autosomes via the supplied index, independently of the gzip scan.
            result = subprocess.run(
                ["bcftools", "query", "-r", REGIONS, "-f", "%CHROM\t%POS\t%REF\t%ALT\n", str(vcf)],
                capture_output=True, timeout=600,
            )
            index_ok = result.returncode == 0 and hashlib.sha256(result.stdout).hexdigest() == digest
            if not index_ok:
                issues.add("indexed_retrieval_failed_or_differs")
            conflicts = 0
            for coordinate, alleles in sites.items():
                if coordinate in cohort_alleles and cohort_alleles[coordinate] != alleles:
                    conflicts += 1
                else:
                    cohort_alleles[coordinate] = alleles
            if conflicts:
                issues.add("cross_sample_ref_alt_conflict")
            if reference_sites is None:
                reference_sites = sites.copy()
            same_sites = sites == reference_sites
            missing_sites = sum(key not in sites for key in reference_sites)
            extra_sites = sum(key not in reference_sites for key in sites)
            rows.append({
                "table_sample_id": sample[id_column], "vcf_sample_id": name,
                "workflow_id": check["workflow_id"], "records": records,
                "expected_records": expected, "indexed_retrieval_matches": index_ok,
                "same_sites_and_alleles_as_first": same_sites,
                "missing_coordinates_vs_first": missing_sites,
                "extra_coordinates_vs_first": extra_sites,
                "allele_conflicts_vs_prior_samples": conflicts,
                "issues": issues, "review": review,
            })
            vcf.unlink()
            index.unlink()
            if number % 25 == 0 or number == EXPECTED_SAMPLES:
                print(f"Checked samples: {number}/{EXPECTED_SAMPLES}", flush=True)

    # Part 4. Reports; preserve all previous outputs
    for indices in sample_name_rows.values():
        if len(indices) > 1:
            for i in indices:
                rows[i]["issues"].add("duplicate_vcf_sample_name")
    for row in rows:
        row["status"] = "FAIL" if row["issues"] else "REVIEW" if row["review"] else "PASS"
        row["issues"] = ";".join(sorted(row["issues"]))
    failures = sum(row["status"] == "FAIL" for row in rows)
    reviews = sum(bool(row["review"]) for row in rows)
    overall = "FAIL" if failures else "REVIEW" if reviews else "PASS"
    summary = [
        f"Overall: {overall}", f"Samples checked: {len(rows)} / {EXPECTED_SAMPLES}",
        f"Samples with failed checks: {failures}", f"Sample ID strings requiring review: {reviews}",
        f"Unique internal sample names: {len(sample_name_rows)}",
        f"Samples matching first sample's sites and alleles: {sum(r['same_sites_and_alleles_as_first'] for r in rows)}",
        f"Samples with allele conflicts: {sum(r['allele_conflicts_vs_prior_samples'] > 0 for r in rows)}",
        "Site-set differences alone do not fail this check. No VCFs are modified.",
        "This does not verify biological identity, reference build or genotype concordance.",
    ]
    REPORT_DIR.mkdir(parents=True, exist_ok=False)
    with (REPORT_DIR / "content_checks.tsv").open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with (REPORT_DIR / "summary.txt").open("x") as handle:
        handle.write("\n".join(summary) + "\n")
    print("\n" + "\n".join(summary))
    print(f"Reports saved in: {REPORT_DIR}")
    return 0 if overall == "PASS" else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except RuntimeError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
    except (OSError, ValueError, EOFError, subprocess.TimeoutExpired) as error:
        print(f"ERROR: {type(error).__name__}; check input files, tools and cloud access.", file=sys.stderr)
        sys.exit(2)
