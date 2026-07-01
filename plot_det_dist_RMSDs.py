#!/usr/bin/env python3

import os
import sys
import subprocess
import re
from collections import defaultdict
from tqdm import tqdm
import matplotlib.pyplot as plt
import plotly.graph_objects as go


def find_datasets(base_dir):
    """
    Find all directories that contain BOTH:
      - refined.expt
      - dials.refine.log
    """
    datasets = []
    for root, dirs, files in os.walk(base_dir):
        if "refined.expt" in files and "dials.refine.log" in files:
            rel = os.path.relpath(root, base_dir)
            datasets.append((rel, root))
    return datasets


def extract_distances(expt_path):
    """
    Run `dials.python.show` and extract detector distances.
    Handles multiple panels by returning a list.
    """
    try:
        result = subprocess.run(
            ["dials.show", expt_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
    except Exception as e:
        print(f"[WARNING] Failed to run dials.show on {expt_path}: {e}")
        return []

    stdout = result.stdout

    # Extract all "distance: XXX" occurrences inside detector blocks
    distances = []
    for match in re.finditer(r"\bdistance:\s*([0-9.+-eE]+)", stdout):
        distances.append(float(match.group(1)))

    distances = [d for d in distances if d > 0]
    return distances

def extract_final_rmsds(log_path):
    """
    Extract the FINAL RMSDs table from refine.log.
    """
    with open(log_path, "r") as f:
        lines = f.readlines()

    # Find all indices where RMSD tables start
    start_indices = [
        i for i, line in enumerate(lines)
        if "RMSDs by experiment" in line
    ]

    if not start_indices:
        return []

    # Take the last occurrence
    start = start_indices[-1]

    # Extract table block (until blank line or non-table content)
    table_lines = []
    for line in lines[start:]:
        if line.strip() == "":
            break
        if "|" in line or "+" in line or "RMSDs by experiment" in line:
            table_lines.append(line)
        elif table_lines:
            break  # stop once table has started and ends

    # Now parse rows
    rmsd_data = []
    row_pattern = re.compile(
        r"\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*([0-9.]+)\s*\|\s*([0-9.]+)\s*\|\s*([0-9.]+)\s*\|"
    )

    for line in table_lines:
        match = row_pattern.search(line)
        if match:
            rmsd_data.append({
                "exp_id": int(match.group(1)),
                "nref": int(match.group(2)),
                "rmsd_x": float(match.group(3)),
                "rmsd_y": float(match.group(4)),
                "rmsd_z": float(match.group(5)),
            })

    return rmsd_data

def assemble_data(base_dir):
    """
    Collect distances + RMSDs across all datasets.
    """
    datasets = find_datasets(base_dir)

    all_data = []

    for label, path in tqdm(datasets):
        expt = os.path.join(path, "refined.expt")
        log = os.path.join(path, "dials.refine.log")

        distances = extract_distances(expt)
        rmsds = extract_final_rmsds(log)

        # Match experiments and panels (best-effort pairing)
        n = max(len(distances), len(rmsds))

        for i in range(n):
            entry = {
                "label": f"{label}_{i}" if n > 1 else label,
                "distance": distances[i] if i < len(distances) else None,
                "rmsd": rmsds[i] if i < len(rmsds) else None,
            }
            all_data.append(entry)

    return all_data


def plot_results(data):
    labels_dist = []
    distances = []

    labels_rx = []
    rmsd_x = []

    labels_ry = []
    rmsd_y = []

    for d in data:
        # Distance plot: only needs distance
        if d["distance"] is not None and d["distance"] > 0:
            labels_dist.append(d["label"])
            distances.append(d["distance"])

        # RMSD plots: only need rmsd
        if d["rmsd"] is not None:
            labels_rx.append(d["label"])
            rmsd_x.append(d["rmsd"]["rmsd_x"])

            labels_ry.append(d["label"])
            rmsd_y.append(d["rmsd"]["rmsd_y"])

    # ---- Distance plot ----
    if distances:
        plt.figure()
        plt.scatter(labels_dist, distances, marker="o")
        plt.xticks(rotation=60, ha="right")
        plt.xlabel("Dataset")
        plt.ylabel("Detector Distance")
        plt.title("Detector Distance by Dataset")
        plt.tight_layout()
        plt.savefig("distance.png")
        plt.close()
    else:
        print("[WARNING] No distance data to plot")

    # ---- RMSD_X plot ----
    if rmsd_x:
        plt.figure()
        plt.scatter(labels_rx, rmsd_x, marker="o")
        plt.xticks(rotation=60, ha="right")
        plt.xlabel("Dataset")
        plt.ylabel("RMSD_X (px)")
        plt.title("RMSD_X by Dataset")
        plt.tight_layout()
        plt.savefig("rmsd_x.png")
        plt.close()
    else:
        print("[WARNING] No RMSD_X data to plot")

    # ---- RMSD_Y plot ----
    if rmsd_y:
        plt.figure()
        plt.scatter(labels_ry, rmsd_y, marker="o")
        plt.xticks(rotation=60, ha="right")
        plt.xlabel("Dataset")
        plt.ylabel("RMSD_Y (px)")
        plt.title("RMSD_Y by Dataset")
        plt.tight_layout()
        plt.savefig("rmsd_y.png")
        plt.close()
    else:
        print("[WARNING] No RMSD_Y data to plot")

    print("[INFO] Plotting complete")


def get_prefix(label):
    return label.split("/")[0]

def get_scan_number(label):
    m = re.search(r"x_(\d+)", label)
    return int(m.group(1)) if m else - 1

def plot_metric_by_prefix(data, metric, ylabel, output_file):
    """
    metric:
        "distance"
        "rmsd_x"
        "rmsd_y"
    """

    grouped = defaultdict(list)

    # ---- group data ----
    for d in data:
        if metric == "distance" and d["distance"] is not None:
            grouped[get_prefix(d["label"])].append(d)

        elif metric in ("rmsd_x", "rmsd_y") and d["rmsd"] is not None:
            grouped[get_prefix(d["label"])].append(d)

    fig = go.Figure()

    # ---- build traces ----
    for prefix, entries in grouped.items():
        # sort by scan number
        entries = sorted(entries, key=lambda d: get_scan_number(d["label"]))

        x = list(range(len(entries)))
        labels = [d["label"] for d in entries]

        if metric == "distance":
            y = [d["distance"] for d in entries]
            hover_value = "Distance"

        elif metric == "rmsd_x":
            y = [d["rmsd"]["rmsd_x"] for d in entries]
            hover_value = "RMSD_X"

        elif metric == "rmsd_y":
            y = [d["rmsd"]["rmsd_y"] for d in entries]
            hover_value = "RMSD_Y"

        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode="markers",
            name=prefix,          # colour by prefix (one trace per prefix)
            text=labels,
            hovertemplate=
                "Dataset: %{text}<br>"
                "Index: %{x}<br>"
                f"{hover_value}: "+"%{y:.3f}<extra></extra>"
        ))

    fig.update_layout(
        title=f"{ylabel} (grouped by prefix)",
        xaxis_title="Dataset index (sorted by scan number)",
        yaxis_title=ylabel,
    )

    fig.write_html(output_file)
    print(f"[INFO] Wrote {output_file}")

def main():
    if len(sys.argv) != 2:
        print("Usage: script.py <base_directory>")
        sys.exit(1)

    base_dir = sys.argv[1]

    data = assemble_data(base_dir)

    print("Collected data:")
    for d in data:
        print(d)

    #plot_results(data)
    plot_metric_by_prefix(data, "distance", "Detector Distance", "distance.html")
    plot_metric_by_prefix(data, "rmsd_x", "RMSD_X (px)", "rmsd_x.html")
    plot_metric_by_prefix(data, "rmsd_y", "RMSD_Y (px)", "rmsd_y.html")


if __name__ == "__main__":
    main()