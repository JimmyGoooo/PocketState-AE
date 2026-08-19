# ColabFold input package

- Reference: `sp|P35408|PE2R4_HUMAN Prostaglandin E2 receptor EP4 subtype OS=Homo sapiens OX=9606 GN=PTGER4 PE=1 SV=1`
- Wild-type sequence length: 488 aa
- Model objective: compare full-length PTGER4 WT with explicitly selected single-site variants.
- Experimental pocket evidence: `experimental_pocket_candidate_sites.csv`, derived from 6 Å heavy-atom contacts to 7UR (5YWY) and P2E (7D7M).

## Files

- `PTGER4_WT.fasta`: validated WT input for ColabFold.
- `PTGER4_WT_and_variant_panel.fasta`: WT plus explicitly selected variants, when present.
- `experimental_pocket_candidate_sites.csv`: candidate locations, not recommended mutations.
- `selected_variants.csv`: user decision sheet. Set `include_in_colabfold` to `YES` only after recording a rationale.

## Generate a variant panel

1. Fill one or more rows in `selected_variants.csv`.
2. Run `python scripts/prepare_colabfold_inputs.py` from the project root.
3. Submit either `PTGER4_WT.fasta` or `PTGER4_variant_panel.fasta` to ColabFold.

ColabFold predictions are sequence-based; preserve 5YWY and 7D7M as the evidence for ligand-bound inactive/active pockets.
