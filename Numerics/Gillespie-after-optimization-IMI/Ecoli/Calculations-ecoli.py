# =============================================================================
# FILE: Numerics/Gillespie-after-optimization-IMI/Ecoli/Calculations-ecoli.py
# PAPER: An information-theoretic perspective on feed-forward loop abundances in transcriptional networks
# AUTHORS: Mintu Nandi, Sudip Chattopadhyay, and Suman K Banik
# CONTACT: mintunandi@ubi.s.u-tokyo.ac.jp; sudip@chem.iiests.ac.in;
#          skbanik@jcbose.ac.in
#
# PURPOSE
#   Combines the Ecoli Gillespie outputs for each motif and gate.
#   For every motif, the script reads the FFL total MI and the matched
#   open-loop joint MI and calculates sMI = MI_FFL - MI_OL.
#   In manuscript notation, sMI = I_int(X;Z).
#
# INPUTS
#   A directory containing all files named
#   total-MI-abund-<motif>ffl-<gate>-num.dat and
#   total-MI-abund-<motif>ol-<gate>-num.dat for c1-c4, i1-i4, and/or.
#   Each input file must contain: <motif label><whitespace><numeric value>.
#
# EXECUTION
#   cd Numerics/Gillespie-after-optimization-IMI/Ecoli
#   python Calculations-ecoli.py --dir .
#
# OUTPUTS
#   sMI-coherent-and-Ecoli-num.dat
#   sMI-incoherent-and-Ecoli-num.dat
#   sMI-coherent-or-Ecoli-num.dat
#   sMI-incoherent-or-Ecoli-num.dat
#
# CODE-TO-MANUSCRIPT NOTATION
#   v1 = I_FFL(X;Z); v2 = I_OL(X,X_tilde;Z)=I_path(X;Z).
#   new_value = v1 - v2 = I_int(X;Z).
#   group_name = coherent or incoherent; gate = and or or.
#
# DEPENDENCIES / NOTES
#   Uses only the Python standard library (pathlib and argparse).
#   Motif labels are checked strictly before subtraction.
# =============================================================================

from pathlib import Path
import argparse

def find_existing_file(base_dir: Path, stem: str) -> Path:
    """
    Tries both with and without .dat, because some systems hide extensions in the file explorer.
    stem example: 'total-MI-abund-c1ffl-and-num'
    """
    candidates = [
        base_dir / stem,
        base_dir / f"{stem}.dat",
    ]
    for p in candidates:
        if p.exists() and p.is_file():
            return p
    raise FileNotFoundError(f"Missing file (tried): {', '.join(str(c) for c in candidates)}")

def read_motif_value(path: Path):
    """
    Expects first non-empty line like:
        C1    0.123456
    Returns (motif_str, value_float)
    """
    with path.open("r") as f:
        for line in f:
            line = line.strip()
            if line:
                parts = line.split()
                if len(parts) < 2:
                    raise ValueError(f"Bad format in {path}: '{line}'")
                return parts[0], float(parts[1])
    raise ValueError(f"Empty file: {path}")

def compute_group(base_dir: Path, motifs, gate: str, group_name: str):
    """
    For each motif:
      total-MI-abund-<motif_lower>ffl-<gate>-num(.dat)
      total-MI-abund-<motif_lower>ol-<gate>-num(.dat)
    Computes: ffl - ol
    Writes:
      sMI-<group_name>-<gate>-Ecoli-num.dat
    """
    out_path = base_dir / f"sMI-{group_name}-{gate}-Ecoli-num.dat"
    out_lines = []

    for motif in motifs:
        motif_lower = motif.lower()  # C1->c1, I1->i1

        ffl_stem = f"total-MI-abund-{motif_lower}ffl-{gate}-num"
        ol_stem  = f"total-MI-abund-{motif_lower}ol-{gate}-num"

        ffl_path = find_existing_file(base_dir, ffl_stem)
        ol_path  = find_existing_file(base_dir, ol_stem)

        m1, v1 = read_motif_value(ffl_path)
        m2, v2 = read_motif_value(ol_path)

        # Strict consistency checks (remove if you prefer)
        if m1.upper() != motif.upper():
            raise ValueError(f"Motif label mismatch: {ffl_path} contains '{m1}', expected '{motif}'")
        if m2.upper() != motif.upper():
            raise ValueError(f"Motif label mismatch: {ol_path} contains '{m2}', expected '{motif}'")

        new_value = v1 - v2
        out_lines.append(f"{motif}\t{new_value:.6f}\n")

    with out_path.open("w") as f:
        f.writelines(out_lines)

    return out_path

def main():
    parser = argparse.ArgumentParser(description="Compute sMI = (ffl MI) - (ol MI) for coherent/incoherent, and/or.")
    parser.add_argument("--dir", type=str, default=".", help="Directory containing the total-MI-abund-* files")
    args = parser.parse_args()

    base_dir = Path(args.dir)

    coherent = ["C1", "C2", "C3", "C4"]
    incoherent = ["I1", "I2", "I3", "I4"]
    gates = ["and", "or"]

    for gate in gates:
        compute_group(base_dir, coherent, gate=gate, group_name="coherent")
        compute_group(base_dir, incoherent, gate=gate, group_name="incoherent")

if __name__ == "__main__":
    main()
