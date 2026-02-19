#!/usr/bin/env python3
import csv
import glob
import os
import re
import shutil
import sys


def extract(text: str, pattern: str) -> str | None:
    match = re.search(pattern, text, re.MULTILINE)
    if not match:
        return None
    return match.group(1).strip()


def load_metrics(glob_pattern: str) -> dict[str, tuple[str, str, str]]:
    metrics: dict[str, tuple[str, str, str]] = {}
    for file_path in sorted(glob.glob(glob_pattern)):
        scenario = os.path.basename(file_path).removesuffix(".time.txt")
        text = open(file_path, encoding="utf-8", errors="replace").read()

        elapsed = extract(
            text,
            r"^\s*Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([^\n]+)$",
        )
        cpu_percent = extract(
            text,
            r"^\s*Percent of CPU this job got:\s*([0-9]+)%\s*$",
        )
        max_rss_kb = extract(
            text,
            r"^\s*Maximum resident set size \(kbytes\):\s*([0-9]+)\s*$",
        )

        if not (elapsed and cpu_percent and max_rss_kb):
            raise RuntimeError(f"failed to parse metrics from {file_path}")

        metrics[scenario] = (elapsed, cpu_percent, max_rss_kb)

    return metrics


def main() -> int:
    python_metrics = load_metrics("tmp/python_parity_timings_parallel8/*.time.txt")
    rust_metrics = load_metrics("tmp/rust_parity_timings_parallel8/*.time.txt")

    scenarios = sorted(set(python_metrics) | set(rust_metrics))
    output = "tmp/parity_resource_metrics_parallel8.tsv"

    with open(output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "scenario",
                "python_elapsed",
                "python_cpu_percent",
                "python_max_rss_kb",
                "rust_elapsed",
                "rust_cpu_percent",
                "rust_max_rss_kb",
            ]
        )
        for scenario in scenarios:
            py = python_metrics.get(scenario, ("", "", ""))
            rs = rust_metrics.get(scenario, ("", "", ""))
            writer.writerow([scenario, *py, *rs])

    shutil.copyfile(output, "tmp/parity_resource_metrics.tsv")
    print(output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
