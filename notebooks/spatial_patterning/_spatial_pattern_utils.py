"""Shared helpers for spatial patterning demos (hex lattice + clone labeling).

Standalone repo layout (run notebooks from this directory)::

    .
    ├── data/cell_anno.tsv          # optional input (see load_or_make_barcodes_by_type)
    ├── output/<experiment>/...     # written by notebooks
    ├── _spatial_pattern_utils.py
    └── 01_vary_shape.ipynb / 02_... / 03_...

Requires: numpy, pandas. Optional: matplotlib for plotting in notebooks.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

TARGET_PER_TYPE = {"normal": 1000, "tumor1": 300, "tumor2": 700}
N_TARGET = sum(TARGET_PER_TYPE.values())
SPOT_SPACING = 1.0
SEED = 12345

# Default paths relative to the demo-repo root
DEFAULT_CELL_ANNO = Path("data") / "cell_anno.tsv"
DEFAULT_OUTPUT_DIR = Path("output")

PALETTE = {"normal": "#7f7f7f", "tumor1": "#2ca02c", "tumor2": "#d62728"}

SHAPE_DISPLAY = {
    "well_separated_clusters": "Separated_clusters",
    "mixed_unstructured": "Mixed",
    "ring": "Ring",
    "stripes": "Stripes",
    "gradient_infiltration": "Intermixing",
    "single_tumor_region": "Single_tumor_region",
}


def resolve_repo_root(start: Optional[Path] = None) -> Path:
    """Directory that contains ``_spatial_pattern_utils.py`` (demo-repo root)."""
    start = Path.cwd().resolve() if start is None else Path(start).resolve()
    for cand in [start, *start.parents]:
        if (cand / "_spatial_pattern_utils.py").exists():
            return cand
    return start


def make_hex_grid(n_spots: int = N_TARGET, spacing: float = SPOT_SPACING) -> pd.DataFrame:
    """Visium-like odd-r hexagonal lattice, centered at the origin."""
    side = int(np.ceil(np.sqrt(n_spots)))
    dy = (np.sqrt(3) / 2.0) * spacing
    rows = []
    spot_i = 0
    for r in range(side):
        x_shift = 0.5 * spacing if (r % 2 == 1) else 0.0
        for c in range(side):
            x = c * spacing + x_shift
            y = r * dy
            rows.append((spot_i, r, c, x, y))
            spot_i += 1
            if spot_i >= n_spots:
                break
        if spot_i >= n_spots:
            break
    df = pd.DataFrame(rows, columns=["spot_id", "row", "col", "x", "y"])
    df["x"] = df["x"] - df["x"].mean()
    df["y"] = df["y"] - df["y"].mean()
    return df


def _select_top_n(
    df: pd.DataFrame,
    score: np.ndarray,
    n: int,
    exclude_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    if exclude_mask is None:
        exclude_mask = np.zeros(len(df), dtype=bool)
    avail = np.where(~exclude_mask)[0]
    if n > len(avail):
        raise ValueError(f"Not enough available spots: need {n}, have {len(avail)}")
    order = avail[np.argsort(score[avail])[::-1]]
    return order[:n]


def assign_exact_counts(
    df: pd.DataFrame, tumor2_idx: np.ndarray, tumor1_idx: np.ndarray
) -> pd.Series:
    labels = pd.Series("normal", index=df.index)
    labels.iloc[tumor2_idx] = "tumor2"
    labels.iloc[tumor1_idx] = "tumor1"
    vc = labels.value_counts().to_dict()
    for ct, n in TARGET_PER_TYPE.items():
        if vc.get(ct, 0) != n:
            raise AssertionError(f"Count mismatch for {ct}: got {vc.get(ct, 0)}, expected {n}")
    return labels


# ---- Shape patterns ----

def pattern_mixed_unstructured(df: pd.DataFrame, seed: int = SEED) -> pd.Series:
    r = np.random.default_rng(seed)
    idx = np.arange(len(df))
    r.shuffle(idx)
    tumor2_idx = idx[: TARGET_PER_TYPE["tumor2"]]
    tumor1_idx = idx[
        TARGET_PER_TYPE["tumor2"] : TARGET_PER_TYPE["tumor2"] + TARGET_PER_TYPE["tumor1"]
    ]
    return assign_exact_counts(df, tumor2_idx=tumor2_idx, tumor1_idx=tumor1_idx)


def pattern_well_separated_clusters(df: pd.DataFrame) -> pd.Series:
    x, y = df["x"].to_numpy(), df["y"].to_numpy()
    c1 = (-12.0 * SPOT_SPACING, 0.0)
    c2 = (12.0 * SPOT_SPACING, 0.0)
    d1 = (x - c1[0]) ** 2 + (y - c1[1]) ** 2
    d2 = (x - c2[0]) ** 2 + (y - c2[1]) ** 2
    tumor2_idx = _select_top_n(df, -d2, TARGET_PER_TYPE["tumor2"])
    mask2 = np.zeros(len(df), dtype=bool)
    mask2[tumor2_idx] = True
    tumor1_idx = _select_top_n(df, -d1, TARGET_PER_TYPE["tumor1"], exclude_mask=mask2)
    return assign_exact_counts(df, tumor2_idx=tumor2_idx, tumor1_idx=tumor1_idx)


def pattern_ring(df: pd.DataFrame) -> pd.Series:
    x, y = df["x"].to_numpy(), df["y"].to_numpy()
    r2 = np.sqrt(x**2 + y**2)
    score_t2 = -np.abs(r2 - 7.0 * SPOT_SPACING)
    score_t1 = -np.abs(r2 - 4.0 * SPOT_SPACING)
    tumor2_idx = _select_top_n(df, score_t2, TARGET_PER_TYPE["tumor2"])
    mask2 = np.zeros(len(df), dtype=bool)
    mask2[tumor2_idx] = True
    tumor1_idx = _select_top_n(df, score_t1, TARGET_PER_TYPE["tumor1"], exclude_mask=mask2)
    return assign_exact_counts(df, tumor2_idx=tumor2_idx, tumor1_idx=tumor1_idx)


def pattern_stripes(df: pd.DataFrame) -> pd.Series:
    y = df["y"].to_numpy()
    y0_t2 = 1.25 * SPOT_SPACING
    y0_t1 = -1.75 * SPOT_SPACING
    y_split = (y0_t1 + y0_t2) / 2.0
    tumor2_idx = _select_top_n(
        df, -np.abs(y - y0_t2), TARGET_PER_TYPE["tumor2"], exclude_mask=(y < y_split)
    )
    exclude_t1 = (y > y_split) | np.isin(np.arange(len(df)), tumor2_idx)
    tumor1_idx = _select_top_n(
        df, -np.abs(y - y0_t1), TARGET_PER_TYPE["tumor1"], exclude_mask=exclude_t1
    )
    return assign_exact_counts(df, tumor2_idx=tumor2_idx, tumor1_idx=tumor1_idx)


def pattern_gradient_infiltration(df: pd.DataFrame, seed: int = SEED) -> pd.Series:
    """Fixed cross-infiltration (15% tumor1 / 10% tumor2) for the Intermixing shape."""
    x, y = df["x"].to_numpy(), df["y"].to_numpy()
    c1 = (-10.0 * SPOT_SPACING, 0.0)
    c2 = (10.0 * SPOT_SPACING, 0.0)
    d1 = (x - c1[0]) ** 2 + (y - c1[1]) ** 2
    d2 = (x - c2[0]) ** 2 + (y - c2[1]) ** 2
    r = np.random.default_rng(seed)

    left_need = TARGET_PER_TYPE["tumor1"] + int(round(0.10 * TARGET_PER_TYPE["tumor2"]))
    right_need = TARGET_PER_TYPE["tumor2"] + int(round(0.15 * TARGET_PER_TYPE["tumor1"]))
    left_pool = _select_top_n(df, -d1, min(len(df), int(round(left_need * 1.25))))
    right_pool = _select_top_n(df, -d2, min(len(df), int(round(right_need * 1.15))))

    overlap = np.intersect1d(left_pool, right_pool)
    if overlap.size:
        closer_to_left = d1[overlap] < d2[overlap]
        left_pool = np.concatenate([np.setdiff1d(left_pool, overlap), overlap[closer_to_left]])
        right_pool = np.concatenate([np.setdiff1d(right_pool, overlap), overlap[~closer_to_left]])

    n_t1_in_right = int(round(0.15 * TARGET_PER_TYPE["tumor1"]))
    n_t2_in_left = int(round(0.10 * TARGET_PER_TYPE["tumor2"]))
    t1_in_right = r.choice(right_pool, size=n_t1_in_right, replace=False)
    t2_in_left = r.choice(left_pool, size=n_t2_in_left, replace=False)
    t2_main = r.choice(
        np.setdiff1d(right_pool, t1_in_right),
        size=TARGET_PER_TYPE["tumor2"] - n_t2_in_left,
        replace=False,
    )
    t1_main = r.choice(
        np.setdiff1d(left_pool, t2_in_left),
        size=TARGET_PER_TYPE["tumor1"] - n_t1_in_right,
        replace=False,
    )
    tumor2_idx = np.concatenate([t2_main, t2_in_left])
    tumor1_idx = np.concatenate([t1_main, t1_in_right])
    if np.intersect1d(tumor1_idx, tumor2_idx).size:
        raise AssertionError("tumor1/tumor2 overlap in gradient_infiltration")
    return assign_exact_counts(df, tumor2_idx=tumor2_idx, tumor1_idx=tumor1_idx)


def pattern_single_tumor_region(df: pd.DataFrame) -> pd.Series:
    x, y = df["x"].to_numpy(), df["y"].to_numpy()
    c2 = (5.0 * SPOT_SPACING, 0.0)
    c1 = (2.0 * SPOT_SPACING, -3.0 * SPOT_SPACING)
    d2 = (x - c2[0]) ** 2 + (y - c2[1]) ** 2
    d1 = (x - c1[0]) ** 2 + (y - c1[1]) ** 2
    tumor2_idx = _select_top_n(df, -d2, TARGET_PER_TYPE["tumor2"])
    mask2 = np.zeros(len(df), dtype=bool)
    mask2[tumor2_idx] = True
    tumor1_idx = _select_top_n(df, -d1, TARGET_PER_TYPE["tumor1"], exclude_mask=mask2)
    return assign_exact_counts(df, tumor2_idx=tumor2_idx, tumor1_idx=tumor1_idx)


SHAPE_FUNCS = {
    "well_separated_clusters": pattern_well_separated_clusters,
    "mixed_unstructured": pattern_mixed_unstructured,
    "ring": pattern_ring,
    "stripes": pattern_stripes,
    "gradient_infiltration": pattern_gradient_infiltration,
    "single_tumor_region": pattern_single_tumor_region,
}


# ---- Distance / mixing ----

def pattern_vary_distance(df: pd.DataFrame, distance: float) -> pd.Series:
    y = df["y"].to_numpy()
    half_gap = distance / 2.0
    y0_t2 = half_gap
    y0_t1 = -half_gap
    y_split = 0.0
    tumor2_idx = _select_top_n(
        df, -np.abs(y - y0_t2), TARGET_PER_TYPE["tumor2"], exclude_mask=(y < y_split)
    )
    exclude_t1 = (y > y_split) | np.isin(np.arange(len(df)), tumor2_idx)
    tumor1_idx = _select_top_n(
        df, -np.abs(y - y0_t1), TARGET_PER_TYPE["tumor1"], exclude_mask=exclude_t1
    )
    return assign_exact_counts(df, tumor2_idx=tumor2_idx, tumor1_idx=tumor1_idx)


def pattern_vary_mixing_rate(
    df: pd.DataFrame, infiltration_rate: float, seed: int = SEED
) -> pd.Series:
    x, y = df["x"].to_numpy(), df["y"].to_numpy()
    c1 = (-10.0 * SPOT_SPACING, 0.0)
    c2 = (10.0 * SPOT_SPACING, 0.0)
    d1 = (x - c1[0]) ** 2 + (y - c1[1]) ** 2
    d2 = (x - c2[0]) ** 2 + (y - c2[1]) ** 2
    rng = np.random.default_rng(seed)

    n_t1_in_right = int(round(infiltration_rate * TARGET_PER_TYPE["tumor1"]))
    n_t2_in_left = int(round(infiltration_rate * TARGET_PER_TYPE["tumor2"]))
    n_t1_main = TARGET_PER_TYPE["tumor1"] - n_t1_in_right
    n_t2_main = TARGET_PER_TYPE["tumor2"] - n_t2_in_left

    left_need = n_t1_main + n_t2_in_left
    right_need = n_t2_main + n_t1_in_right
    left_pool = _select_top_n(df, -d1, min(len(df), int(round(left_need * 1.3))))
    right_pool = _select_top_n(df, -d2, min(len(df), int(round(right_need * 1.3))))

    overlap = np.intersect1d(left_pool, right_pool)
    if overlap.size:
        closer_to_left = d1[overlap] < d2[overlap]
        left_pool = np.concatenate([np.setdiff1d(left_pool, overlap), overlap[closer_to_left]])
        right_pool = np.concatenate([np.setdiff1d(right_pool, overlap), overlap[~closer_to_left]])

    t1_in_right = (
        rng.choice(right_pool, size=n_t1_in_right, replace=False)
        if n_t1_in_right
        else np.array([], dtype=int)
    )
    t2_in_left = (
        rng.choice(left_pool, size=n_t2_in_left, replace=False)
        if n_t2_in_left
        else np.array([], dtype=int)
    )
    t2_main = rng.choice(np.setdiff1d(right_pool, t1_in_right), size=n_t2_main, replace=False)
    t1_main = rng.choice(np.setdiff1d(left_pool, t2_in_left), size=n_t1_main, replace=False)

    tumor2_idx = np.concatenate([t2_main, t2_in_left])
    tumor1_idx = np.concatenate([t1_main, t1_in_right])
    if np.intersect1d(tumor1_idx, tumor2_idx).size:
        raise AssertionError("tumor1/tumor2 overlap")
    return assign_exact_counts(df, tumor2_idx=tumor2_idx, tumor1_idx=tumor1_idx)


def distance_schedule(grid: pd.DataFrame, n: int = 5) -> np.ndarray:
    y_min, y_max = grid["y"].min(), grid["y"].max()
    d_max = y_max - y_min - 2.0
    return np.linspace(0, d_max, n)


def mixing_rate_schedule(n: int = 6, rate_max: float = 0.4) -> np.ndarray:
    return np.linspace(0, rate_max, n)


# ---- Barcodes + I/O ----
#
# INPUT  (optional)
#   data/cell_anno.tsv
#     - format: header-free TSV, two columns: barcode <TAB> clone_label
#     - clone_label in {normal, tumor1, tumor2}
#     - need >= 1000 normal, >= 300 tumor1, >= 700 tumor2 rows
#     - if missing, synthetic barcodes normal_0000 / tumor1_0000 / ... are used
#
# OUTPUT (per condition directory under output/<experiment>/<condition>/)
#   tissue_positions_list.csv
#     - Space Ranger–style, NO header
#     - columns: barcode, in_tissue, x, y, pixel_row, pixel_col
#   spot_anno_pattern.tsv
#     - TSV WITH header; index = barcode; column spot_anno in {normal,tumor1,tumor2}
#   barcodes.tsv.gz
#     - one barcode per line (gzipped), same order as positions file


def load_or_make_barcodes_by_type(
    cell_anno_path: Optional[Path] = None,
    seed: int = SEED,
) -> Dict[str, np.ndarray]:
    """Load barcodes from ``data/cell_anno.tsv`` or synthesize placeholders.

    Parameters
    ----------
    cell_anno_path :
        Path to header-free TSV with columns ``barcode``, ``clone_label``
        (``normal`` / ``tumor1`` / ``tumor2``). If ``None`` or missing, synthetic
        IDs are generated.
    """
    rng = np.random.default_rng(seed)
    if cell_anno_path is not None and Path(cell_anno_path).exists():
        df = pd.read_csv(cell_anno_path, sep="\t", header=None, names=["barcode", "spot_anno"])
        out = {}
        for ct, n in TARGET_PER_TYPE.items():
            pool = df.loc[df["spot_anno"] == ct, "barcode"].to_numpy()
            if len(pool) < n:
                raise ValueError(f"Need >= {n} barcodes for {ct}; found {len(pool)}")
            idx = rng.permutation(len(pool))[:n]
            out[ct] = pool[idx]
        return out

    out = {}
    for ct, n in TARGET_PER_TYPE.items():
        out[ct] = np.array([f"{ct}_{i:04d}" for i in range(n)], dtype=object)
    return out


def map_barcodes_to_labels(
    labels: pd.Series,
    barcodes_by_type: Dict[str, np.ndarray],
    seed: int,
) -> List[str]:
    rng = np.random.default_rng(seed)
    assigned = np.empty(len(labels), dtype=object)
    for ct in TARGET_PER_TYPE:
        pos = np.where(labels.values == ct)[0]
        bc = barcodes_by_type[ct].copy()
        rng.shuffle(bc)
        if len(bc) != len(pos):
            raise ValueError(f"{ct}: {len(bc)} barcodes vs {len(pos)} spots")
        assigned[pos] = bc
    return assigned.tolist()


def write_pattern_dir(
    out_dir: Path,
    df_grid: pd.DataFrame,
    labels: pd.Series,
    barcodes: Sequence[str],
) -> Path:
    """Write Space Ranger–style spatial outputs for one condition.

    Files written under ``out_dir``:

    - ``tissue_positions_list.csv`` — no header; ``barcode,in_tissue,x,y,pixel_row,pixel_col``
    - ``spot_anno_pattern.tsv`` — barcode index + ``spot_anno`` column
    - ``barcodes.tsv.gz`` — one barcode per line
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if len(barcodes) != len(df_grid):
        raise ValueError(f"barcodes length {len(barcodes)} != n_spots {len(df_grid)}")

    pos = pd.DataFrame(
        {
            "barcode": list(barcodes),
            "in_tissue": 1,
            "x": df_grid["x"].values,
            "y": df_grid["y"].values,
            "pixel_row": np.round(df_grid["x"].values).astype(int),
            "pixel_col": np.round(df_grid["y"].values).astype(int),
        }
    )
    pos.to_csv(out_dir / "tissue_positions_list.csv", header=False, index=False)

    anno = pd.DataFrame(
        {"spot_anno": labels.values}, index=pd.Index(list(barcodes), name="barcode")
    )
    anno.to_csv(out_dir / "spot_anno_pattern.tsv", sep="\t")

    with gzip.open(out_dir / "barcodes.tsv.gz", "wt") as f:
        for b in barcodes:
            f.write(f"{b}\n")
    return out_dir


def plot_spatial(
    ax,
    df_grid: pd.DataFrame,
    labels: pd.Series,
    title: str = "",
    spot_size: float = 8,
):
    for ct, color in PALETTE.items():
        mask = labels.values == ct
        if not mask.any():
            continue
        ax.scatter(
            df_grid.loc[mask, "x"],
            df_grid.loc[mask, "y"],
            s=spot_size,
            c=color,
            label=ct,
            linewidths=0,
        )
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    if title:
        ax.set_title(title, fontsize=10)
