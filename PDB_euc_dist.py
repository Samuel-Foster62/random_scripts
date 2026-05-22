#!/usr/bin/env python3
"""
import sys
import math
import numpy as np


# -----------------------------
# CONFIG
# -----------------------------
IGNORE_HYDROGENS = True
STRICT_MATCH = False  # If True → fail on missing atoms


# -----------------------------
# PDB PARSER
# -----------------------------
def parse_pdb_atoms(pdb_file):
    atoms = {}

    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                if IGNORE_HYDROGENS and atom_name.startswith("H"):
                    continue

                altloc = line[16].strip()
                resname = line[17:20].strip()
                chain = line[21].strip()
                resseq = int(line[22:26])
                icode = line[26].strip()

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                occ = float(line[54:60])

                key = (chain, resseq, icode, resname, atom_name)

                # Handle altloc: keep highest occupancy
                if key not in atoms or occ > atoms[key]["occ"]:
                    atoms[key] = {
                        "coords": np.array([x, y, z]),
                        "line": line,
                        "occ": occ,
                        "altloc": altloc
                    }

    return atoms


# -----------------------------
# DISTANCE
# -----------------------------
def compute_distances(coords1, coords2):
    diffs = coords1 - coords2
    return np.linalg.norm(diffs, axis=1)


# -----------------------------
# KABSCH ALIGNMENT
# -----------------------------
def kabsch(P, Q):
    
    #Align Q onto P
    #Returns rotated Q and RMSD
    
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)

    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q

    C = np.dot(Q_centered.T, P_centered)
    V, S, Wt = np.linalg.svd(C)

    d = np.sign(np.linalg.det(np.dot(V, Wt)))
    D = np.diag([1, 1, d])

    U = np.dot(V, np.dot(D, Wt))

    Q_rot = np.dot(Q_centered, U)

    Q_aligned = Q_rot + centroid_P

    rmsd = np.sqrt(np.mean(np.sum((P - Q_aligned) ** 2, axis=1)))

    return Q_aligned, rmsd


# -----------------------------
# MAIN
# -----------------------------
def main(pdb1, pdb2, output_pdb):

    atoms1 = parse_pdb_atoms(pdb1)
    atoms2 = parse_pdb_atoms(pdb2)

    keys1 = set(atoms1.keys())
    keys2 = set(atoms2.keys())

    common_keys = sorted(keys1 & keys2)

    missing_in_1 = keys2 - keys1
    missing_in_2 = keys1 - keys2

    # Report missing atoms
    if missing_in_1:
        print(f"⚠ Missing in {pdb1}: {len(missing_in_1)} atoms")
    if missing_in_2:
        print(f"⚠ Missing in {pdb2}: {len(missing_in_2)} atoms")

    if STRICT_MATCH and (missing_in_1 or missing_in_2):
        raise ValueError("Strict mode enabled and mismatches found.")

    if not common_keys:
        raise ValueError("No common atoms found!")

    # Build coordinate arrays
    coords1 = np.array([atoms1[k]["coords"] for k in common_keys])
    coords2 = np.array([atoms2[k]["coords"] for k in common_keys])

    # RMSD before alignment
    pre_rmsd = np.sqrt(np.mean(np.sum((coords1 - coords2) ** 2, axis=1)))
    print(f"Initial RMSD: {pre_rmsd:.3f} Å")

    # Align
    coords2_aligned, post_rmsd = kabsch(coords1, coords2)
    print(f"Aligned RMSD: {post_rmsd:.3f} Å")

    # Compute distances after alignment
    distances = compute_distances(coords1, coords2_aligned)

    # Map distances back to keys
    dist_dict = dict(zip(common_keys, distances))

    # Write output
    with open(pdb1) as fin, open(output_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                if IGNORE_HYDROGENS and atom_name.startswith("H"):
                    fout.write(line)
                    continue

                resname = line[17:20].strip()
                chain = line[21].strip()
                resseq = int(line[22:26])
                icode = line[26].strip()

                key = (chain, resseq, icode, resname, atom_name)

                if key in dist_dict:
                    d = dist_dict[key]
                else:
                    d = 0.0  # fallback for missing

                new_b = f"{d:6.2f}"
                new_line = line[:60] + new_b + line[66:]
                fout.write(new_line)
            else:
                fout.write(line)

    print(f"✅ Output written to {output_pdb}")


# -----------------------------
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python pdb_dist.py file1.pdb file2.pdb output.pdb")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3])

"""
#!/usr/bin/env python3

import sys
import numpy as np
import math
import matplotlib.pyplot as plt


IGNORE_HYDROGENS = True


# -----------------------------
# PDB PARSING
# -----------------------------
def parse_pdb(pdb_file):
    residues = {}
    atoms = {}

    with open(pdb_file) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip()
                if IGNORE_HYDROGENS and atom_name.startswith("H"):
                    continue

                resname = line[17:20].strip()
                chain = line[21].strip()
                resseq = int(line[22:26])
                icode = line[26].strip()

                key_res = (chain, resseq, icode)
                key_atom = (chain, resseq, icode, resname, atom_name)

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                residues[key_res] = resname
                atoms[key_atom] = np.array([x, y, z])

    return residues, atoms


# -----------------------------
# NEEDLEMAN-WUNSCH
# -----------------------------
def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
    n, m = len(seq1), len(seq2)

    score = np.zeros((n+1, m+1))
    traceback = np.zeros((n+1, m+1), dtype=int)

    for i in range(1, n+1):
        score[i,0] = i * gap
    for j in range(1, m+1):
        score[0,j] = j * gap

    for i in range(1, n+1):
        for j in range(1, m+1):
            match_score = score[i-1, j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete = score[i-1, j] + gap
            insert = score[i, j-1] + gap

            score[i,j] = max(match_score, delete, insert)

            if score[i,j] == match_score:
                traceback[i,j] = 0
            elif score[i,j] == delete:
                traceback[i,j] = 1
            else:
                traceback[i,j] = 2

    # traceback
    align1, align2 = [], []
    i, j = n, m

    while i > 0 or j > 0:
        if i > 0 and j > 0 and traceback[i,j] == 0:
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and traceback[i,j] == 1:
            align1.append(seq1[i-1])
            align2.append("-")
            i -= 1
        else:
            align1.append("-")
            align2.append(seq2[j-1])
            j -= 1

    return align1[::-1], align2[::-1]


# -----------------------------
# KABSCH
# -----------------------------
def kabsch(P, Q):
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)

    P -= centroid_P
    Q -= centroid_Q

    H = Q.T @ P
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1,1,d]) @ Vt

    Q_rot = Q @ R
    Q_aligned = Q_rot + centroid_P

    rmsd = np.sqrt(np.mean(np.sum((P + centroid_P - Q_aligned)**2, axis=1)))

    return Q_aligned, rmsd


# -----------------------------
def main(pdb1, pdb2, output_pdb):
    res1, atoms1 = parse_pdb(pdb1)
    res2, atoms2 = parse_pdb(pdb2)

    keys1 = sorted(res1.keys())
    keys2 = sorted(res2.keys())

    seq1 = [res1[k] for k in keys1]
    seq2 = [res2[k] for k in keys2]

    aln1, aln2 = needleman_wunsch(seq1, seq2)

    # map residues
    mapping = []
    i = j = 0

    for a, b in zip(aln1, aln2):
        if a != "-" and b != "-":
            mapping.append((keys1[i], keys2[j]))
            i += 1
            j += 1
        elif a != "-":
            i += 1
        else:
            j += 1

    # build atom pairs
    coords1 = []
    coords2 = []
    atom_keys = []

    for r1, r2 in mapping:
        for atom_name in ["CA"]:  # use CA for alignment
            key1 = (*r1, res1[r1], atom_name)
            key2 = (*r2, res2[r2], atom_name)

            if key1 in atoms1 and key2 in atoms2:
                coords1.append(atoms1[key1])
                coords2.append(atoms2[key2])

    coords1 = np.array(coords1)
    coords2 = np.array(coords2)

    # superpose
    coords2_aligned, rmsd = kabsch(coords1.copy(), coords2.copy())

    print(f"Aligned RMSD: {rmsd:.3f} Å")

    # full atom distances
    distances = []
    residue_index = []

    for idx, (r1, r2) in enumerate(mapping):
        for atom_name in ["CA"]:
            key1 = (*r1, res1[r1], atom_name)
            key2 = (*r2, res2[r2], atom_name)

            if key1 in atoms1 and key2 in atoms2:
                d = np.linalg.norm(atoms1[key1] - coords2_aligned[idx])
                distances.append(d)
                residue_index.append(r1[1])

    # -----------------------------
    # PLOTS
    # -----------------------------
    plt.figure()
    plt.plot(residue_index, distances)
    plt.xlabel("Residue Number")
    plt.ylabel("Distance (Å)")
    plt.title("Distance per Residue")
    plt.savefig("distance_plot.png", dpi=300)

    plt.figure()
    plt.hist(distances, bins=50)
    plt.xlabel("Distance (Å)")
    plt.ylabel("Frequency")
    plt.title("Distance Distribution")
    plt.savefig("distance_histogram.png", dpi=300)

    print("✅ Plots saved:")
    print(" - distance_plot.png")
    print(" - distance_histogram.png")


# -----------------------------
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py file1.pdb file2.pdb output.pdb")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2], sys.argv[3])