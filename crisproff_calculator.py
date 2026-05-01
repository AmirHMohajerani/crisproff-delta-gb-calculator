#!/usr/bin/env python3
"""
CRISPRoff Binding Energy Calculator (Delta G_B)
================================================
Thermodynamic model for gRNA:DNA binding energy based on:
  - Sugimoto et al. (1995) Biochemistry 34:11211      — RNA:DNA NN parameters
  - Watkins & SantaLucia (2011) NAR 39:1894           — RNA:DNA mismatch NN parameters
  - SantaLucia (1998) PNAS 95:1460                    — DNA:DNA NN parameters (Delta G_O)
  - Alkan et al. (2018) Genome Biology 19:177         — Cas9 positional weights Gamma[i]

Formula:
  Delta G_B = delta_PAM x (Delta G_H + 3.1 - Delta G_U - Delta G_O)

  Delta G_H = sum over 19 steps of Gamma[i] x NN[i]
    - matched step  : Sugimoto RNA:DNA NN value
    - mismatched step: Watkins mismatch NN value (type-specific)
  Delta G_O = sum of SantaLucia DNA:DNA NN values for target sequence
  Delta G_U = gRNA self-folding MFE (set to 0 if not calculated)
  delta_PAM : 1.0 (NGG), 0.9 (NAG), 0.8 (NGA)

Usage:
  python crisproff_calculator.py

Author: built with AI assistance (Perplexity) as part of MSc thesis work
        on allele-specific CRISPRoff epigenetic silencing — A.H. Mohajerani Araghi
"""

# ── Parameter tables ──────────────────────────────────────────────────────────

# Sugimoto 1995: RNA:DNA nearest-neighbour Delta G°37 (kcal/mol), matched steps
# Key = gRNA dinucleotide (5'→3'), assumes Watson-Crick complementarity on DNA
SUGIMOTO_NN = {
    'AA': -1.0,  'AC': -2.1,  'AG': -1.8,  'AT': -0.9,
    'CA': -0.9,  'CC': -2.1,  'CG': -1.7,  'CT': -0.9,
    'GA': -1.3,  'GC': -2.7,  'GG': -2.9,  'GT': -1.1,
    'TA': -0.6,  'TC': -1.6,  'TG': -1.4,  'TT': -0.2,
}

# Watkins & SantaLucia 2011: RNA:DNA mismatch NN Delta G°37 (kcal/mol)
# Key = gRNA_base + complement(DNA_coding_base)  →  mismatch type in RNA:DNA hybrid
# Self-mismatches: Table 1 dimer averages (experimentally measured)
# Cross-mismatches: estimates based on SantaLucia 1998 DNA:DNA data
MISMATCH_NN = {
    'AA': +0.90,   # rA:dA  — Watkins 2011 Table 1 avg
    'CC': +1.20,   # rC:dC  — Watkins 2011 Table 1 avg
    'GG': +0.08,   # rG:dG  — Watkins 2011 Table 1 avg (most tolerated same-base)
    'TT': +0.59,   # rU:dT  — Watkins 2011 Table 1 avg
    'GT': +0.30,   # rG:dT wobble — most tolerated mismatch (Xiang 2019)
    'GA': +0.60,   # rG:dA  — estimate
    'TG': +0.70,   # rU:dG  — estimate
    'TC': +0.90,   # rU:dC  — estimate
    'AG': +1.20,   # rA:dG  — estimate
    'AC': +1.30,   # rA:dC  — estimate
    'CA': +1.20,   # rC:dA  — estimate
    'CT': +1.10,   # rC:dT  — estimate
}

# SantaLucia 1998: DNA:DNA nearest-neighbour Delta G°37 (kcal/mol)
# Used for Delta G_O (energy cost of opening the target DNA duplex)
SANTALUCIA_NN = {
    'AA': -1.00,  'AT': -0.88,  'TA': -0.58,  'CA': -1.45,
    'GT': -1.44,  'CT': -1.28,  'GA': -1.30,  'CG': -2.17,
    'GC': -2.24,  'GG': -1.84,  'TT': -1.00,  'TG': -1.44,
    'AC': -1.45,  'AG': -1.28,  'TC': -1.30,  'CC': -1.84,
}

# Alkan 2018: Cas9 positional weights Gamma[i] for 19 dinucleotide steps
# REVERSED so index 0 = PAM-distal (low weight), index 18 = PAM-proximal (high weight)
# Matches 5'→3' gRNA entry with PAM at the 3' end
ALKAN_WEIGHTS = list(reversed([
    2.25, 2.04, 2.15, 1.94, 2.08, 1.93, 2.02, 1.72, 1.88,
    1.98, 1.51, 1.39, 1.00, 1.46, 1.38, 2.13, 1.90, 1.96, 1.80
]))

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
PAM_FACTOR = {'NGG': 1.0, 'NAG': 0.9, 'NGA': 0.8}


# ── Core calculation ──────────────────────────────────────────────────────────

def calculate_dgb(grna: str, dna: str, pam: str = 'NGG', dg_u: float = 0.0,
                  verbose: bool = False) -> dict:
    """
    Calculate Delta G_B for a gRNA:DNA pair.

    Parameters
    ----------
    grna    : 20-nt gRNA sequence as DNA (T not U), 5'→3'
    dna     : 20-nt target DNA protospacer (same strand as gRNA, non-template), 5'→3'
    pam     : PAM type — 'NGG' (default), 'NAG', or 'NGA'
    dg_u    : gRNA self-folding MFE from RNAfold (kcal/mol); use 0 if unknown
    verbose : print step-by-step breakdown

    Returns
    -------
    dict with keys: dg_h, dg_o, dg_u, dg_b, score, mismatches, mismatch_positions
    """
    grna = grna.upper().strip()
    dna  = dna.upper().strip()

    if len(grna) != 20 or len(dna) != 20:
        raise ValueError("Both gRNA and DNA must be exactly 20 nucleotides.")

    mismatch_positions = [i+1 for i in range(20) if grna[i] != dna[i]]

    if verbose:
        print(f"\ngRNA : 5'-{grna}-3'")
        print(f"DNA  : 5'-{dna}-3'")
        print(f"Mismatches at positions (1-based): {mismatch_positions or 'none'}")
        print(f"\n{'Step':>4}  {'gRNA_din':>8}  {'Match':>5}  {'NN_val':>8}  "
              f"{'Gamma':>6}  {'Contrib':>9}")
        print("-" * 50)

    dg_h = 0.0
    for i in range(19):
        is_mismatch = (grna[i] != dna[i])
        if is_mismatch:
            mm_key = grna[i] + COMPLEMENT.get(dna[i], 'N')
            nn_val = MISMATCH_NN.get(mm_key, +1.0)
            flag   = f"MM({mm_key})"
        else:
            nn_val = SUGIMOTO_NN.get(grna[i] + grna[i+1], 0.0)
            flag   = "✓"
        contrib = ALKAN_WEIGHTS[i] * nn_val
        dg_h   += contrib
        if verbose:
            print(f"{i+1:>4}  {grna[i]+grna[i+1]:>8}  {flag:<10}  "
                  f"{nn_val:>8.2f}  {ALKAN_WEIGHTS[i]:>6.2f}  {contrib:>9.4f}")

    dg_h_init = dg_h + 3.1
    dg_o      = sum(SANTALUCIA_NN.get(dna[i] + dna[i+1], 0.0) for i in range(19))
    delta_pam = PAM_FACTOR.get(pam.upper(), 1.0)
    dg_b      = delta_pam * (dg_h_init - dg_u - dg_o)

    if verbose:
        print(f"\n  Delta G_H                 = {dg_h:.4f} kcal/mol")
        print(f"  + initiation (+3.1)       = {dg_h_init:.4f} kcal/mol")
        print(f"  Delta G_O (DNA duplex)    = {dg_o:.4f} kcal/mol")
        print(f"  Delta G_U (gRNA fold)     = {dg_u:.4f} kcal/mol")
        print(f"  delta_PAM ({pam})          = {delta_pam}")
        print(f"  Delta G_B                 = {dg_b:.4f} kcal/mol")
        print(f"  CRISPRoff score (-DG_B)   = {-dg_b:.4f}")

    return {
        'dg_h':               round(dg_h, 4),
        'dg_h_init':          round(dg_h_init, 4),
        'dg_o':               round(dg_o, 4),
        'dg_u':               dg_u,
        'dg_b':               round(dg_b, 4),
        'score':              round(-dg_b, 4),
        'mismatches':         len(mismatch_positions),
        'mismatch_positions': mismatch_positions,
        'binds':              dg_b <= -10.0,
    }


def compare_alleles(grna: str, dna_mut: str, dna_wt: str,
                    pam: str = 'NGG', dg_u: float = 0.0) -> None:
    """
    Compare gRNA binding to a mutant and wild-type allele.
    Prints a summary table and ΔΔG_B.
    """
    print("=" * 60)
    print("  ALLELE SPECIFICITY ANALYSIS")
    print("=" * 60)

    r_mut = calculate_dgb(grna, dna_mut, pam, dg_u, verbose=False)
    r_wt  = calculate_dgb(grna, dna_wt,  pam, dg_u, verbose=False)
    ddg   = r_wt['dg_b'] - r_mut['dg_b']

    print(f"  gRNA          : 5'-{grna}-3'")
    print(f"  MUT target    : 5'-{dna_mut}-3'")
    print(f"  WT target     : 5'-{dna_wt}-3'")
    print(f"  PAM           : {pam}")
    print()
    print(f"  {'':30} {'MUT':>10} {'WT':>10}")
    print(f"  {'Mismatches':30} {r_mut['mismatches']:>10} {r_wt['mismatches']:>10}")
    print(f"  {'Delta G_B (kcal/mol)':30} {r_mut['dg_b']:>10.3f} {r_wt['dg_b']:>10.3f}")
    print(f"  {'CRISPRoff score':30} {r_mut['score']:>10.3f} {r_wt['score']:>10.3f}")
    print(f"  {'Binds? (threshold -10)':30} {'YES' if r_mut['binds'] else 'NO':>10} "
          f"{'YES' if r_wt['binds'] else 'NO':>10}")
    print()
    print(f"  ΔΔG_B (WT - MUT) = {ddg:+.3f} kcal/mol")
    if ddg > 2.0 and not r_wt['binds']:
        print("  ✓ ALLELE-SPECIFIC: gRNA discriminates MUT from WT")
    elif ddg > 2.0 and r_wt['binds']:
        print("  ⚠ Thermodynamic difference present but WT still binds.")
        print("    Consider adding a deliberate PAM-proximal mismatch.")
    else:
        print("  ✗ Insufficient discrimination. Redesign gRNA.")
    print("=" * 60)


# ── Example run ──────────────────────────────────────────────────────────────

if __name__ == "__main__":

    # Your thesis sequences
    GRNA    = "TCTTAAAGATAATGATAGCC"   # designed to match MUT allele
    DNA_MUT = "TCTTAAAGATAATGATAGCC"   # mutant allele (perfect match)
    DNA_WT  = "TGTTAAAGCTAATGATAGCC"   # wild-type allele (1 PAM-distal mismatch)

    print("\n── Single calculation (verbose) ──")
    result = calculate_dgb(GRNA, DNA_MUT, pam="NGG", verbose=True)

    print("\n── Allele specificity comparison ──")
    compare_alleles(GRNA, DNA_MUT, DNA_WT, pam="NGG")
