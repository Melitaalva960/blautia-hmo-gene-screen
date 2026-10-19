# Marker gene panel

Save the protein sequences of your marker panel here as `hmo_panel.faa` (FASTA).
Use short, unique IDs in the headers, because they become the column names of the
presence/absence matrix, for example:

```
>Blon_2335_LNBP
MSKH...
>BbrUCC2003_RhGLnbp112
MTEA...
```

The thesis used:

- **HMO panel:** 64 genes across 19 pathways, taken from *Bifidobacterium breve*,
  *B. longum*, *B. bifidum* and *Roseburia hominis*, classified with CAZy.
- **FOS/GOS panel:** 10 sequences: five GH32 enzymes (FOS) from *Blautia faecicola*,
  *B. obeum* and *B. massiliensis*; four beta-galactosidases and one ABC transporter
  solute-binding protein (GOS) from *B. breve*.

Add `panel_metadata.tsv` (gene ID, pathway, function, organism, UniProt/NCBI accession)
so others can see where each sequence came from.
