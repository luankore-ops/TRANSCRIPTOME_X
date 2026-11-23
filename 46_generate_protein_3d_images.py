
"""
TRANSCRIPTOME-X
SCRIPT 46 - PROTEIN 3D IMAGE GENERATION

Generates:
- High-resolution PNG images
- Interactive HTML 3D visualizations
- Structural audit CSV

Sources:
- Original AlphaFold models
- Recovered AlphaFold models
- Experimental PDB structures

Original structural files are never modified.
"""

from __future__ import annotations

import csv
import logging
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

import plotly.graph_objects as go

from Bio.PDB import MMCIFParser, PDBParser


# ============================================================
# CONFIGURATION
# ============================================================

PROJECT_DIR = Path(
    r"C:\Users\luank\Dropbox\PC\Desktop\PYTHON\TRANSCRIPTOME_X"
)

BASE_DIR = (
    PROJECT_DIR
    / "results"
    / "gene_protein_illustrations"
)

OUTPUT_DIR = BASE_DIR / "final_3d_images"

PNG_DIR = OUTPUT_DIR / "png"
HTML_DIR = OUTPUT_DIR / "html"
AUDIT_DIR = OUTPUT_DIR / "audit"

# Search all these directories recursively.
SEARCH_DIRECTORIES = [
    BASE_DIR / "alphafold_3d",
    BASE_DIR / "cache",
    BASE_DIR / "structure_recovery" / "alphafold",
    BASE_DIR / "structure_recovery" / "experimental_pdb",
]

GENES = [
    ("THBS1", "P07996"),
    ("CCN3", "P48745"),
    ("NRP1", "O14786"),
    ("TGFBR3", "Q03167"),
    ("VEGFC", "P49767"),
    ("PIK3CA", "P42336"),
    ("FGF1", "P05230"),
    ("RELN", "P78509"),
    ("MET", "P08581"),
    ("PAK2", "Q13177"),
    ("COL4A2", "P08572"),
    ("FBLN1", "P23142"),
    ("CCN2", "P29279"),
    ("COL5A1", "P20908"),
    ("PLOD2", "O00469"),
    ("LUM", "P51884"),
    ("MMP1", "P03956"),
    ("MATN2", "O00339"),
    ("ADAMTS6", "Q9UKP5"),
    ("HTRA1", "Q92743"),
    ("FRAS1", "Q86XX4"),
    ("TAGLN", "Q01995"),
    ("TIMP4", "Q99727"),
]

SUPPORTED_EXTENSIONS = {".cif", ".mmcif", ".pdb", ".ent"}

PNG_DPI = 300
FIGURE_SIZE = (12, 10)

MAX_PLOTLY_POINTS = 15000
MAX_MATPLOTLIB_POINTS = 25000

# ============================================================
# DIRECTORIES AND LOGGING
# ============================================================

for directory in [
    PNG_DIR,
    HTML_DIR,
    AUDIT_DIR,
]:
    directory.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[
        logging.FileHandler(
            AUDIT_DIR / "script46.log",
            encoding="utf-8",
        ),
        logging.StreamHandler(sys.stdout),
    ],
)

logger = logging.getLogger("TRANSCRIPTOME-X")


# ============================================================
# STRUCTURE DISCOVERY
# ============================================================

def identify_source(path: Path) -> str:

    normalized = str(path).lower()

    if "experimental_pdb" in normalized:
        return "EXPERIMENTAL_PDB"

    if "structure_recovery" in normalized:
        return "RECOVERED_ALPHAFOLD"

    return "ORIGINAL_ALPHAFOLD"


def find_structure_files(
    gene: str,
    accession: str,
) -> list[Path]:

    candidates = []

    identifiers = [
        gene.lower(),
        accession.lower(),
    ]

    for directory in SEARCH_DIRECTORIES:

        if not directory.exists():
            continue

        for path in directory.rglob("*"):

            if not path.is_file():
                continue

            if path.suffix.lower() not in SUPPORTED_EXTENSIONS:
                continue

            if path.stat().st_size < 100:
                continue

            filename = path.name.lower()

            if any(identifier in filename for identifier in identifiers):
                candidates.append(path)

    # Also inspect original project structure files whose names
    # may contain only the gene symbol.
    unique = list(dict.fromkeys(candidates))

    def priority(path: Path):

        source = identify_source(path)

        order = {
            "ORIGINAL_ALPHAFOLD": 0,
            "RECOVERED_ALPHAFOLD": 1,
            "EXPERIMENTAL_PDB": 2,
        }

        return (
            order.get(source, 9),
            -path.stat().st_size,
        )

    return sorted(unique, key=priority)


# ============================================================
# STRUCTURE PARSING
# ============================================================

def parse_structure(path: Path, gene: str):

    suffix = path.suffix.lower()

    with warnings.catch_warnings():

        warnings.simplefilter("ignore")

        if suffix in {".cif", ".mmcif"}:

            parser = MMCIFParser(
                QUIET=True,
                auth_chains=True,
            )

        else:

            parser = PDBParser(
                QUIET=True,
                PERMISSIVE=True,
            )

        structure = parser.get_structure(
            gene,
            str(path),
        )

    chain_data = []

    for model in structure:

        for chain in model:

            coordinates = []
            b_factors = []
            residue_numbers = []

            for residue in chain:

                if residue.id[0] != " ":
                    continue

                if not residue.has_id("CA"):
                    continue

                atom = residue["CA"]

                coordinates.append(
                    np.asarray(atom.coord, dtype=float)
                )

                b_factors.append(
                    float(atom.get_bfactor())
                )

                residue_numbers.append(
                    int(residue.id[1])
                )

            if len(coordinates) < 2:
                continue

            chain_data.append({
                "chain": str(chain.id),
                "coordinates": np.asarray(coordinates),
                "b_factors": np.asarray(b_factors),
                "residue_numbers": np.asarray(residue_numbers),
            })

    if not chain_data:
        raise ValueError(
            "No valid C-alpha coordinates found."
        )

    return chain_data


# ============================================================
# DATA PREPARATION
# ============================================================

def prepare_plot_data(chain_data):

    all_coordinates = np.vstack([
        item["coordinates"]
        for item in chain_data
    ])

    center = np.mean(all_coordinates, axis=0)

    centered = []

    for item in chain_data:

        coordinates = item["coordinates"] - center

        centered.append({
            **item,
            "coordinates": coordinates,
        })

    return centered


def calculate_structure_metrics(chain_data):

    total_residues = sum(
        len(item["coordinates"])
        for item in chain_data
    )

    total_chains = len(chain_data)

    all_bfactors = np.concatenate([
        item["b_factors"]
        for item in chain_data
    ])

    return {
        "residues": total_residues,
        "chains": total_chains,
        "mean_bfactor": float(np.mean(all_bfactors)),
        "median_bfactor": float(np.median(all_bfactors)),
        "min_bfactor": float(np.min(all_bfactors)),
        "max_bfactor": float(np.max(all_bfactors)),
    }


# ============================================================
# DOWNSAMPLING
# ============================================================

def sample_indices(length: int, maximum: int):

    if length <= maximum:
        return np.arange(length)

    return np.unique(
        np.linspace(
            0,
            length - 1,
            maximum,
            dtype=int,
        )
    )


# ============================================================
# MATPLOTLIB PNG GENERATION
# ============================================================

def generate_png(
    gene: str,
    accession: str,
    chain_data,
    source: str,
    output_path: Path,
):

    fig = plt.figure(figsize=FIGURE_SIZE)

    ax = fig.add_subplot(
        111,
        projection="3d",
    )

    all_values = np.concatenate([
        item["b_factors"]
        for item in chain_data
    ])

    color_label = (
        "pLDDT confidence"
        if source != "EXPERIMENTAL_PDB"
        else "Experimental B-factor"
    )

    vmin = float(np.min(all_values))
    vmax = float(np.max(all_values))

    if vmax <= vmin:
        vmax = vmin + 1.0

    scatter = None

    for item in chain_data:

        coords = item["coordinates"]
        values = item["b_factors"]

        indices = sample_indices(
            len(coords),
            MAX_MATPLOTLIB_POINTS,
        )

        sampled = coords[indices]
        sampled_values = values[indices]

        # Backbone trace.
        ax.plot(
            sampled[:, 0],
            sampled[:, 1],
            sampled[:, 2],
            linewidth=0.65,
            alpha=0.45,
        )

        scatter = ax.scatter(
            sampled[:, 0],
            sampled[:, 1],
            sampled[:, 2],
            c=sampled_values,
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
            s=5,
            alpha=0.9,
        )

    if scatter is not None:

        colorbar = fig.colorbar(
            scatter,
            ax=ax,
            shrink=0.65,
            pad=0.08,
        )

        colorbar.set_label(color_label)

    ax.set_title(
        f"{gene} | {accession}\n"
        f"3D C-alpha backbone | {source}",
        fontsize=13,
        pad=18,
    )

    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_zlabel("Z (Å)")

    ax.set_box_aspect((1, 1, 1))

    ax.view_init(
        elev=25,
        azim=45,
    )

    fig.tight_layout()

    fig.savefig(
        output_path,
        dpi=PNG_DPI,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(fig)


# ============================================================
# PLOTLY HTML GENERATION
# ============================================================

def generate_html(
    gene: str,
    accession: str,
    chain_data,
    source: str,
    output_path: Path,
):

    fig = go.Figure()

    color_label = (
        "pLDDT"
        if source != "EXPERIMENTAL_PDB"
        else "B-factor"
    )

    for item in chain_data:

        coords = item["coordinates"]
        values = item["b_factors"]

        indices = sample_indices(
            len(coords),
            MAX_PLOTLY_POINTS,
        )

        sampled = coords[indices]
        sampled_values = values[indices]

        chain = item["chain"]

        fig.add_trace(
            go.Scatter3d(
                x=sampled[:, 0],
                y=sampled[:, 1],
                z=sampled[:, 2],
                mode="lines+markers",
                name=f"Chain {chain}",
                line={
                    "width": 3,
                    "color": "rgba(120,120,120,0.45)",
                },
                marker={
                    "size": 2.5,
                    "color": sampled_values,
                    "colorscale": "Viridis",
                    "cmin": float(np.min(values)),
                    "cmax": max(
                        float(np.max(values)),
                        float(np.min(values)) + 1.0,
                    ),
                    "colorbar": {
                        "title": color_label,
                    },
                    "showscale": True,
                },
                customdata=np.column_stack([
                    item["residue_numbers"][indices],
                    sampled_values,
                ]),
                hovertemplate=(
                    "Chain: " + chain +
                    "<br>Residue: %{customdata[0]}" +
                    "<br>" + color_label + ": %{customdata[1]:.2f}" +
                    "<extra></extra>"
                ),
            )
        )

    fig.update_layout(
        title={
            "text": (
                f"{gene} ({accession})"
                f"<br><sup>3D C-alpha backbone | {source}</sup>"
            ),
            "x": 0.5,
        },
        scene={
            "xaxis_title": "X (Å)",
            "yaxis_title": "Y (Å)",
            "zaxis_title": "Z (Å)",
            "aspectmode": "data",
        },
        template="plotly_white",
        height=850,
        width=1100,
        legend={
            "orientation": "h",
            "y": -0.05,
        },
    )

    fig.write_html(
        str(output_path),
        include_plotlyjs=True,
        full_html=True,
    )


# ============================================================
# MAIN
# ============================================================

def main():

    audit_records = []

    logger.info("=" * 70)
    logger.info("TRANSCRIPTOME-X - SCRIPT 46")
    logger.info("PROTEIN 3D IMAGE GENERATION")
    logger.info("Total proteins: %s", len(GENES))
    logger.info("=" * 70)

    for index, (gene, accession) in enumerate(
        GENES,
        start=1,
    ):

        logger.info(
            "\n[%s/%s] %s",
            index,
            len(GENES),
            gene,
        )

        record = {
            "gene": gene,
            "uniprot_accession": accession,
            "source": "",
            "structure_file": "",
            "chains": 0,
            "ca_residues": 0,
            "mean_bfactor": "",
            "png_status": "NOT_GENERATED",
            "html_status": "NOT_GENERATED",
            "status": "NOT_FOUND",
            "error": "",
        }

        candidates = find_structure_files(
            gene,
            accession,
        )

        if not candidates:

            logger.warning(
                "%s: no structural file found.",
                gene,
            )

            record["status"] = "NO_STRUCTURE_AVAILABLE"
            audit_records.append(record)

            continue

        selected_path = candidates[0]
        source = identify_source(selected_path)

        record["source"] = source
        record["structure_file"] = str(selected_path)

        try:

            chain_data = parse_structure(
                selected_path,
                gene,
            )

            chain_data = prepare_plot_data(chain_data)

            metrics = calculate_structure_metrics(
                chain_data
            )

            record["chains"] = metrics["chains"]
            record["ca_residues"] = metrics["residues"]
            record["mean_bfactor"] = round(
                metrics["mean_bfactor"],
                3,
            )

            png_path = PNG_DIR / f"{gene}.png"
            html_path = HTML_DIR / f"{gene}.html"

            # Generate PNG.
            generate_png(
                gene,
                accession,
                chain_data,
                source,
                png_path,
            )

            record["png_status"] = (
                "SUCCESS"
                if png_path.exists() and png_path.stat().st_size > 1000
                else "FAILED"
            )

            # Generate interactive HTML.
            generate_html(
                gene,
                accession,
                chain_data,
                source,
                html_path,
            )

            record["html_status"] = (
                "SUCCESS"
                if html_path.exists() and html_path.stat().st_size > 1000
                else "FAILED"
            )

            if (
                record["png_status"] == "SUCCESS"
                and record["html_status"] == "SUCCESS"
            ):
                record["status"] = "SUCCESS"
            else:
                record["status"] = "PARTIAL"

            logger.info(
                "%s: %s C-alpha residues | %s chains",
                gene,
                metrics["residues"],
                metrics["chains"],
            )

            logger.info(
                "%s: PNG=%s | HTML=%s",
                gene,
                record["png_status"],
                record["html_status"],
            )

        except Exception as exc:

            logger.exception(
                "%s: structure visualization failed.",
                gene,
            )

            record["status"] = "ERROR"
            record["error"] = str(exc)

        audit_records.append(record)

    # ========================================================
    # SAVE AUDIT
    # ========================================================

    audit_path = AUDIT_DIR / "protein_3d_image_audit.csv"

    with audit_path.open(
        "w",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=list(audit_records[0].keys()),
        )

        writer.writeheader()
        writer.writerows(audit_records)

    # ========================================================
    # FINAL SUMMARY
    # ========================================================

    png_success = sum(
        item["png_status"] == "SUCCESS"
        for item in audit_records
    )

    html_success = sum(
        item["html_status"] == "SUCCESS"
        for item in audit_records
    )

    structure_found = sum(
        item["status"] in {"SUCCESS", "PARTIAL"}
        for item in audit_records
    )

    print("\n" + "=" * 70)
    print("SCRIPT 46 - COMPLETED")
    print("=" * 70)

    print(f"Total proteins: {len(GENES)}")
    print(f"Structures found: {structure_found}/{len(GENES)}")
    print(f"PNG images generated: {png_success}/{len(GENES)}")
    print(f"Interactive HTML generated: {html_success}/{len(GENES)}")

    print("\nRESULTS BY PROTEIN")

    for item in audit_records:

        print(
            f"{item['gene']:10s} | "
            f"{item['status']:24s} | "
            f"PNG: {item['png_status']:12s} | "
            f"HTML: {item['html_status']}"
        )

    print("\nOUTPUT:")
    print(f"PNG:   {PNG_DIR}")
    print(f"HTML:  {HTML_DIR}")
    print(f"AUDIT: {audit_path}")

    print("=" * 70)


if __name__ == "__main__":
    main()
