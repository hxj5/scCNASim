# Spatial patterning demos

Standalone notebooks that generate Visium-like spatial architectures used in the
manuscript spatial patterning experiments (hexagonal lattice; exact clone counts
1000 normal / 300 tumor1 / 700 tumor2).

**Run all notebooks from this repository root** (the directory that contains
`_spatial_pattern_utils.py`).

| Notebook | Task | Output directory |
|----------|------|------------------|
| `01_vary_shape.ipynb` | Six clonal shapes | `output/vary_shape/` |
| `02_vary_distance.ipynb` | Five interclonal distances | `output/vary_distance/` |
| `03_vary_mixing_rate.ipynb` | Six mixing rates | `output/vary_mixing_rate/` |

Shared code: `_spatial_pattern_utils.py`.

```bash
pip install numpy pandas matplotlib
jupyter nbconvert --to notebook --execute 01_vary_shape.ipynb
jupyter nbconvert --to notebook --execute 02_vary_distance.ipynb
jupyter nbconvert --to notebook --execute 03_vary_mixing_rate.ipynb
```

---

## Directory layout

```
.
├── README.md
├── _spatial_pattern_utils.py
├── 01_vary_shape.ipynb
├── 02_vary_distance.ipynb
├── 03_vary_mixing_rate.ipynb
├── data/
│   ├── README.md
│   └── cell_anno.tsv          # optional input (example included)
└── output/                    # created when notebooks run (gitignored)
    ├── vary_shape/<condition>/
    ├── vary_distance/<condition>/
    └── vary_mixing_rate/<condition>/
```

---

## Input

| Path | Required? | Description |
|------|-----------|-------------|
| `data/cell_anno.tsv` | No | Header-free TSV: `barcode` then `clone_label` |

**Format** (`data/cell_anno.tsv`):

```
AAACCTGAGTTAAGTG-1	normal
AAACCTGGTAGCAAAT-1	normal
...
TCGCGAGCAGATCCAT-1	tumor1
...
```

- Columns (no header): `barcode`, `clone_label`
- `clone_label` ∈ `{normal, tumor1, tumor2}`
- Need at least **1000** normal, **300** tumor1, and **700** tumor2 rows (extras are randomly downsampled)
- If the file is missing, notebooks synthesize barcodes (`normal_0000`, `tumor1_0000`, …)

See [`data/README.md`](data/README.md).

---

## Output

Each condition directory under `output/<experiment>/<condition>/` contains:

| File | Type | Description |
|------|------|-------------|
| `tissue_positions_list.csv` | CSV, **no header** | Space Ranger–style positions: `barcode, in_tissue, x, y, pixel_row, pixel_col` |
| `spot_anno_pattern.tsv` | TSV **with header** | Index = barcode; column `spot_anno` ∈ `{normal, tumor1, tumor2}` |
| `barcodes.tsv.gz` | gzip text | One barcode per line (same order as the positions file) |

Overview figures (`*_overview.png`) are also written under each experiment’s `output/` folder.

These spatial files can be paired with a shared expression matrix / BAM for CalicoST or other spatial tools.
