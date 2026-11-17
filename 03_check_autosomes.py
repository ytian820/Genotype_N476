#!/usr/bin/env python3
"""Check N476 autosomal screening logs and output object availability.

Phenotype: none. Unit: one single-sample workflow in the specified submission.
Inputs: GCS object metadata and stdout; no VCF or index contents are downloaded.
Outputs: workflow_checks.tsv and summary.txt, created only when run by the user.
Rules: 476 workflows; nonnegative integer counts; input equals exclusions plus
retained; retained > 0; completion message present; VCF and index both nonempty.
Scope: checks submission artifacts, not sample-table identity, VCF contents,
index validity, or biological QC. Different retained counts are allowed.
"""

import csv
import json
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


# Part 1. Fixed run configuration
SUBMISSION_ID = "f707a3f7-c2b6-4f36-b8e4-ea477c0681e7"
BUCKET = "gs://fc-secure-6c44e649-ddc0-4b8b-a14d-3ae8f9538b09"
ROOT = f"{BUCKET}/submissions/{SUBMISSION_ID}/ScreenAutosomes/"
EXPECTED_WORKFLOWS = 476
WORKERS = 4
REPORT_DIR = Path("/home/jupyter/CHIP_CAD_N476/result/03_check_autosomes")
COUNT_FIELDS = [
    "input", "excluded_non_autosomal", "excluded_FILTER", "excluded_non_SNP",
    "excluded_duplicate_coordinate_records", "retained",
]
ARTIFACTS = {
    "stdout": "stdout",
    "vcf": "autosomal.screened.vcf.gz",
    "index": "autosomal.screened.vcf.gz.tbi",
}


# Part 2. Read-only cloud access and per-workflow checks
def gcloud(*args):
    result = subprocess.run(
        ["gcloud", "storage", *args], capture_output=True, text=True, timeout=300
    )
    if result.returncode:
        raise RuntimeError(
            "gcloud read failed; check the Terra account, bucket access and network."
        )
    return result.stdout


def list_objects():
    """List sizes once; --long emits bytes, creation time, and GCS URI."""
    objects = {}
    listing = gcloud("ls", "--long", ROOT + "**")
    for line in listing.splitlines():
        if not line.strip() or line.lstrip().startswith("TOTAL:"):
            continue
        match = re.fullmatch(r"\s*(\d+)\s+\S+\s+(gs://.+?)\s*", line)
        if not match or not match[2].startswith(ROOT):
            raise ValueError("Unexpected object listing format; inspection stopped.")
        objects[match[2]] = int(match[1])
    if not objects:
        raise ValueError("No objects found in the specified workflow directory.")
    return objects


def check_workflow(workflow_id, objects):
    prefix = ROOT + workflow_id + "/call-ScreenAutosomalVcf/"
    # Cromwell's first attempt is directly under call-...; retries use attempt-N.
    attempts = set()
    for url in objects:
        if url.startswith(prefix):
            relative = url[len(prefix):]
            if "/" not in relative:
                attempts.add(1)
            else:
                match = re.match(r"attempt-(\d+)/[^/]+$", relative)
                if match:
                    attempts.add(int(match[1]))
    attempt = max(attempts, default=1)
    task_dir = prefix if attempt == 1 else prefix + f"attempt-{attempt}/"
    row = {"workflow_id": workflow_id, "attempt": attempt}
    issues = []
    for label, filename in ARTIFACTS.items():
        url = task_dir + filename
        size = objects.get(url)
        row[label + "_uri"] = url
        row[label + "_bytes"] = "" if size is None else size
        if size is None:
            issues.append("missing_" + label)
        elif size == 0:
            issues.append("empty_" + label)

    if row["stdout_bytes"] not in ("", 0):
        try:
            log = gcloud("cat", row["stdout_uri"])
            lines = [line.removeprefix("Screening_counts:").strip()
                     for line in log.splitlines()
                     if line.startswith("Screening_counts:")]
            if "Autosomal_screening: complete" not in log.splitlines():
                issues.append("missing_completion_message")
            if len(lines) != 1:
                issues.append("expected_one_counts_line")
            else:
                counts = json.loads(lines[0])
                if not isinstance(counts, dict) or any(
                    type(counts.get(key)) is not int or counts[key] < 0
                    for key in COUNT_FIELDS
                ):
                    issues.append("invalid_counts")
                else:
                    row.update({key: counts[key] for key in COUNT_FIELDS})
                    if counts["input"] != sum(counts[key] for key in COUNT_FIELDS[1:]):
                        issues.append("counts_do_not_sum_to_input")
                    if counts["retained"] == 0:
                        issues.append("zero_retained")
        except (RuntimeError, subprocess.TimeoutExpired):
            issues.append("stdout_read_failed")
        except json.JSONDecodeError:
            issues.append("invalid_counts_json")
    row["status"] = "FAIL" if issues else "PASS"
    row["issues"] = ";".join(issues)
    return row


# Part 3. Preflight and batch checking
def main():
    print(f"Submission: {SUBMISSION_ID}", flush=True)
    print(f"Expected workflows: {EXPECTED_WORKFLOWS}; readers: {WORKERS}")
    print(f"Report directory: {REPORT_DIR}")
    print("Checks: stdout counts and nonempty output objects; no VCF content scan.")
    if sys.version_info < (3, 9):
        raise RuntimeError("Python 3.9 or newer is required.")
    if not shutil.which("gcloud"):
        raise RuntimeError("gcloud is required; run this in Terra Terminal.")
    if REPORT_DIR.exists():
        raise RuntimeError("Report directory already exists; nothing overwritten.")

    objects = list_objects()
    workflow_ids = sorted({url[len(ROOT):].split("/")[0] for url in objects})
    # Reject unexpected layouts rather than silently ignoring their objects.
    if any(not re.fullmatch(r"[0-9a-fA-F-]{36}", wid) for wid in workflow_ids):
        raise ValueError("Unexpected workflow directory names; inspection stopped.")
    print(f"Workflow directories found: {len(workflow_ids)}", flush=True)
    rows = []
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        for row in pool.map(lambda wid: check_workflow(wid, objects), workflow_ids):
            rows.append(row)
            if len(rows) % 50 == 0 or len(rows) == len(workflow_ids):
                print(f"Checked workflows: {len(rows)}/{len(workflow_ids)}", flush=True)

    # Part 4. Reports and exit status (no overwrite)
    failed = sum(row["status"] == "FAIL" for row in rows)
    passed = len(rows) - failed
    complete = len(rows) == EXPECTED_WORKFLOWS and failed == 0
    summary = [
        f"Submission: {SUBMISSION_ID}",
        f"Overall: {'PASS' if complete else 'FAIL'}",
        f"Workflow directories: {len(rows)} / {EXPECTED_WORKFLOWS}",
        f"Passed workflows: {passed}",
        f"Failed workflows: {failed}",
    ]
    for label in ARTIFACTS:
        count = sum(isinstance(row[label + "_bytes"], int)
                    and row[label + "_bytes"] > 0 for row in rows)
        summary.append(f"Nonempty {label} objects: {count} / {EXPECTED_WORKFLOWS}")
    valid = [row for row in rows if all(key in row for key in COUNT_FIELDS)]
    summary.append(f"Parsed count records: {len(valid)} / {EXPECTED_WORKFLOWS}")
    if valid:
        summary.append("\nMetric\tMinimum\tMaximum")
        for key in COUNT_FIELDS:
            values = [row[key] for row in valid]
            summary.append(f"{key}\t{min(values)}\t{max(values)}")
    summary.append(
        "\nScope: latest task attempt per workflow; no sample-table identity, "
        "VCF content or index validity checks."
    )

    REPORT_DIR.mkdir(parents=True, exist_ok=False)
    columns = ["workflow_id", "attempt", "status", "issues"] + COUNT_FIELDS
    for label in ARTIFACTS:
        columns.extend([label + "_bytes", label + "_uri"])
    with (REPORT_DIR / "workflow_checks.tsv").open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with (REPORT_DIR / "summary.txt").open("x") as handle:
        handle.write("\n".join(summary) + "\n")
    print("\n" + "\n".join(summary))
    print(f"\nReports saved in: {REPORT_DIR}")
    return 0 if complete else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (RuntimeError, ValueError, OSError, subprocess.TimeoutExpired) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
