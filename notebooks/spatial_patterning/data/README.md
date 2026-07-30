# Input data

## `cell_anno.tsv` (optional)

Header-free tab-separated file mapping barcodes to clone labels.

| Column | Type | Values |
|--------|------|--------|
| 1. `barcode` | string | Unique spot / cell ID (e.g. 10x-style barcode) |
| 2. `clone_label` | string | `normal`, `tumor1`, or `tumor2` |

**Example:**

```
AAACCTGAGTTAAGTG-1	normal
AAACCTGGTCACCTAA-1	tumor1
TCGCGAGCAGATCCAT-1	tumor2
```

**Requirements for demos:**

- At least 1000 rows with `normal`
- At least 300 rows with `tumor1`
- At least 700 rows with `tumor2`

Notebooks randomly subset to exactly those counts (seed `12345`).

If this file is absent, demos still run using synthetic barcodes
(`normal_0000` … `tumor2_0699`).

The included example is adapted from the stCNASim HCC-3–based simulation
cell annotation used in the manuscript spatial patterning experiments.
