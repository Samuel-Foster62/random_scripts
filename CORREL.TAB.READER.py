import sys
import argparse
import io
import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def check_fwf(correl_file):
    with open(correl_file) as f:
        lines = [line.rstrip("\n") for line in f.readlines()[:20]]

    # Get positions of spaces in each line
    space_positions = [set(i for i, c in enumerate(line) if c == " ") for line in lines]

    # Compare consistency
    common_positions = set.intersection(*space_positions)

    print("Common space positions:", sorted(common_positions))

    print("White space loaded")
    df_ws = pd.read_csv(correl_file, sep=r"\s+", header=None)
    print(df_ws.to_markdown())
    
    print("fixed width loaded")
    df_fwf = pd.read_fwf(correl_file, header=None)
    print(df_fwf.to_markdown())

def print_correl_plot(correl_files, columns, output, labels=None):
        print(f"Processing file(s): ")

        all_d = []
        for i, fp in enumerate(correl_files):
            print(f"\t{fp}")
            df = pd.read_fwf(fp, header=None, names=columns)
            all_d.append(df["DMIN"])

            if labels:
                plt.plot(df["DMIN"].apply(angstrom_axis), df["CORREL"], label=labels[i], marker='o')
            else:
                plt.plot(df["DMIN"].apply(angstrom_axis), df["CORREL"], label=fp, marker='o')

        all_d = pd.concat(all_d).drop_duplicates().sort_values(ascending=False)

        step = max(len(all_d) // 5, 1)
        tick_d = pd.concat([
             all_d.iloc[[0]],
             all_d.iloc[::step],
             all_d.iloc[[-1]]
        ]).drop_duplicates()

        tick_positions = tick_d.apply(angstrom_axis)
        tick_labels = [f"{d:.2f}" for d in tick_d]

        plt.xticks(tick_positions, tick_labels)        
        plt.xlabel(r"D$_{min}$ (Å)")
        plt.ylabel("Correlation Coefficient")
        plt.ylim(0, 1)
        plt.legend()
        plt.savefig(output)

def angstrom_axis(n):
    #take a DMIN in angstroms and convert to plot appropriately 
    # (low number higher resolution asymptotically)
    return 1 / (n**2)

def review(correl_files):
    print("Reviewing files")
    print(f"received {len(correl_files)} files")
    print('=' * os.get_terminal_size().columns)
    
    kept_files = []
    labels = []
    print("Type 'y' to KEEP the file, or press Enter to DROP it:\n")
    for fp in correl_files:
        user_choice = input(f"\tKeep '{fp}'? (y/n): ").strip().lower()
        if user_choice == 'y':
            inp = str(input(f"\tPlease enter a new label for dataset '{fp}' in the final plot (press Enter to KEEP the default): ").strip() or f"{fp}")
            if not inp:
                inp = fp
            kept_files.append(fp)
            labels.append(inp)
    print("\nYour final file list: ")
    print(f"\tFile\t| Label")
    for i, kf in enumerate(kept_files):
        print(f"\t{kf}\t {labels[i]}")
    return kept_files, labels

def main():
    p = argparse.ArgumentParser("File reader and plotter for the output of sftools CORREL command")

    p.add_argument('-f', '--files', required=True, nargs='+', help='Space-separated list of files')    #input files
    p.add_argument('-c', "--check", action="store_true", help='check file format of input files')  #reflections post reindexing after phasing
    p.add_argument('-m', "--manual", action="store_true", help='manually review inputted files, select desired files and rename datasets when plotting')
    p.add_argument('-o', "--output", type=str, default="CORREL_plots.png", help='set output plot filename (default: CORREL_plots.png)')
    args = p.parse_args()

    columns = ["NSHELL", "DMIN", "DMAX", "RFACT", "CORREL", "F1", "F2", "NREF"]
    
    if args.files:
        if args.check:
            check_fwf(args.files[0])
        if args.manual:
            files, labels = review(args.files)
            print_correl_plot(files, columns, args.output, labels)
        else:
            print_correl_plot(args.files, columns, args.output)

if __name__ == "__main__":
    main()