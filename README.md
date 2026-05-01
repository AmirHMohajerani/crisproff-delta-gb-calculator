# crisproff-delta-gb-calculator

Thermodynamic calculator for CRISPRoff gRNA:DNA binding energy (\Delta G_B) and allele-specific discrimination (\Delta\Delta G_B).

## Model sources
- Sugimoto et al. 1995 — RNA:DNA nearest-neighbour parameters
- Watkins & SantaLucia 2011 — RNA:DNA mismatch parameters
- SantaLucia 1998 — DNA:DNA nearest-neighbour parameters
- Alkan et al. 2018 — Cas9 positional weights

## Files
- `crisproff_calculator.py` — Python implementation
- `crisproff_calculator_FIXED_FINAL_v8.xlsx` — Excel implementation

## Use
This tool compares gRNA binding against mutant and wild-type alleles to support allele-specific CRISPRoff design.
