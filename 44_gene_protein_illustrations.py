
# ============================================================
# TRANSCRIPTOME-X
# SCRIPT 44 - GENE AND PROTEIN ILLUSTRATION GENERATOR
#
# Generates:
#   1. Conceptual biological illustrations (PNG)
#   2. AlphaFold predicted protein structures (PNG)
#   3. Metadata and download audit (CSV)
#
# Input:
#   23 robust candidate genes identified by TRANSCRIPTOME-X
#
# Output:
#   results/gene_protein_illustrations/
#
# IMPORTANT:
#   - Schematic illustrations are conceptual, not exact
#     representations of protein domains or molecular topology.
#   - 3D structures are computational AlphaFold predictions.
#   - Missing structures are reported, never fabricated.
# ============================================================

from __future__ import annotations

import csv
import io
import logging
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, FancyArrowPatch, Ellipse
import numpy as np
import requests

from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
import plotly.graph_objects as go


# ============================================================
# 1. PROJECT CONFIGURATION
# ============================================================

PROJECT_ROOT = Path(
    r"C:\Users\luank\Dropbox\PC\Desktop\PYTHON\TRANSCRIPTOME_X"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "results"
    / "gene_protein_illustrations"
)

SCHEMATIC_DIR = OUTPUT_ROOT / "schematic"
STRUCTURE_DIR = OUTPUT_ROOT / "alphafold_3d"
CACHE_DIR = OUTPUT_ROOT / "cache"
AUDIT_DIR = OUTPUT_ROOT / "audit"

for directory in [
    SCHEMATIC_DIR,
    STRUCTURE_DIR,
    CACHE_DIR,
    AUDIT_DIR,
]:
    directory.mkdir(parents=True, exist_ok=True)


# ============================================================
# 2. LOGGING
# ============================================================

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(
            AUDIT_DIR / "script_44.log",
            encoding="utf-8",
        ),
    ],
)

logger = logging.getLogger("TRANSCRIPTOME-X")


# ============================================================
# 3. ROBUST GENE PANEL
# ============================================================

@dataclass(frozen=True)
class Gene:
    symbol: str
    category: str
    role: str
    context: str


GENES = [
    Gene(
        "THBS1",
        "extracellular_matrix",
        "Extracellular matrix organization",
        "Extracellular matrix and cell interactions",
    ),
    Gene(
        "CCN3",
        "extracellular_matrix",
        "Cell communication and extracellular signalling",
        "Extracellular signalling environment",
    ),
    Gene(
        "NRP1",
        "membrane_receptor",
        "Cell-surface receptor and signalling co-receptor",
        "Cell membrane and signalling",
    ),
    Gene(
        "TGFBR3",
        "membrane_receptor",
        "TGF-beta receptor-associated signalling",
        "Cell membrane and growth factor signalling",
    ),
    Gene(
        "VEGFC",
        "secreted_signalling",
        "Vascular endothelial growth factor signalling",
        "Extracellular signalling",
    ),
    Gene(
        "PIK3CA",
        "intracellular_signalling",
        "PI3K signalling pathway",
        "Intracellular signal transduction",
    ),
    Gene(
        "FGF1",
        "secreted_signalling",
        "Fibroblast growth factor signalling",
        "Extracellular growth factor signalling",
    ),
    Gene(
        "RELN",
        "extracellular_matrix",
        "Extracellular signalling and cell interactions",
        "Extracellular environment",
    ),
    Gene(
        "MET",
        "membrane_receptor",
        "Receptor tyrosine kinase signalling",
        "Cell membrane and intracellular signalling",
    ),
    Gene(
        "PAK2",
        "intracellular_signalling",
        "Intracellular signalling and cytoskeletal regulation",
        "Intracellular signalling",
    ),
    Gene(
        "COL4A2",
        "extracellular_matrix",
        "Collagen type IV and basement membrane organization",
        "Basement membrane",
    ),
    Gene(
        "FBLN1",
        "extracellular_matrix",
        "Extracellular matrix organization",
        "Extracellular matrix",
    ),
    Gene(
        "CCN2",
        "extracellular_matrix",
        "Extracellular signalling and matrix-associated processes",
        "Extracellular environment",
    ),
    Gene(
        "COL5A1",
        "extracellular_matrix",
        "Collagen type V and fibril organization",
        "Extracellular matrix",
    ),
    Gene(
        "PLOD2",
        "intracellular_signalling",
        "Collagen modification and maturation",
        "Intracellular processing of collagen",
    ),
    Gene(
        "LUM",
        "extracellular_matrix",
        "Extracellular matrix organization",
        "Collagen-rich extracellular matrix",
    ),
    Gene(
        "MMP1",
        "extracellular_matrix",
        "Extracellular matrix degradation",
        "Extracellular matrix remodelling",
    ),
    Gene(
        "MATN2",
        "extracellular_matrix",
        "Extracellular matrix organization",
        "Extracellular matrix",
    ),
    Gene(
        "ADAMTS6",
        "extracellular_matrix",
        "Extracellular matrix-associated protease",
        "Extracellular matrix",
    ),
    Gene(
        "HTRA1",
        "extracellular_matrix",
        "Protease-associated extracellular processes",
        "Extracellular environment",
    ),
    Gene(
        "FRAS1",
        "extracellular_matrix",
        "Basement membrane organization",
        "Basement membrane",
    ),
    Gene(
        "TAGLN",
        "cytoskeleton",
        "Actin-associated cytoskeletal organization",
        "Cellular cytoskeleton",
    ),
    Gene(
        "TIMP4",
        "extracellular_matrix",
        "Regulation of matrix metalloproteinases",
        "Extracellular matrix remodelling",
    ),
]


# ============================================================
# 4. SESSION CONFIGURATION
# ============================================================

SESSION = requests.Session()

SESSION.headers.update({
    "User-Agent": "TRANSCRIPTOME-X/1.0 scientific research script"
})

REQUEST_TIMEOUT = 40
REQUEST_DELAY = 0.5

DPI = 300

FIGURE_SIZE = (10, 7)

# Maximum number of C-alpha atoms displayed in 3D plots.
MAX_STRUCTURE_POINTS = 12000


# ============================================================
# 5. UTILITIES
# ============================================================

def safe_filename(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def save_figure(fig, path: Path) -> None:
    fig.savefig(
        path,
        dpi=DPI,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close(fig)


def request_json(url: str) -> Optional[dict]:
    try:
        response = SESSION.get(
            url,
            timeout=REQUEST_TIMEOUT,
        )

        if response.status_code == 200:
            return response.json()

        logger.warning(
            "HTTP %s while requesting %s",
            response.status_code,
            url,
        )

    except requests.RequestException as exc:
        logger.warning("Request error: %s", exc)

    return None


# ============================================================
# 6. UNIPROT ACCESSION RETRIEVAL
# ============================================================

def get_uniprot_accession(gene_symbol: str) -> Optional[str]:
    """
    Retrieve a reviewed human UniProt accession.

    Search is restricted to Homo sapiens.
    """

    cache_file = CACHE_DIR / f"{gene_symbol}_uniprot.txt"

    if cache_file.exists():
        accession = cache_file.read_text(
            encoding="utf-8"
        ).strip()

        if accession:
            return accession

    query = (
        f"(gene_exact:{gene_symbol}) "
        "AND (organism_id:9606)"
    )

    url = "https://rest.uniprot.org/uniprotkb/search"

    params = {
        "query": query,
        "format": "json",
        "size": 20,
        "fields": "accession,gene_names,reviewed",
    }

    try:
        response = SESSION.get(
            url,
            params=params,
            timeout=REQUEST_TIMEOUT,
        )

        response.raise_for_status()

        data = response.json()

        results = data.get("results", [])

        if not results:
            logger.warning(
                "No UniProt accession found for %s",
                gene_symbol,
            )
            return None

        # Prefer reviewed Swiss-Prot entries.
        selected = None

        for result in results:
            if result.get("entryType") == "UniProtKB reviewed (Swiss-Prot)":
                selected = result
                break

        if selected is None:
            selected = results[0]

        accession = selected.get("primaryAccession")

        if accession:
            cache_file.write_text(
                accession,
                encoding="utf-8",
            )

            return accession

    except requests.RequestException as exc:
        logger.warning(
            "UniProt retrieval failed for %s: %s",
            gene_symbol,
            exc,
        )

    return None


# ============================================================
# 7. ALPHAFOLD STRUCTURE RETRIEVAL
# ============================================================

def get_alphafold_metadata(
    accession: str,
) -> Optional[dict]:

    url = (
        "https://alphafold.ebi.ac.uk/api/prediction/"
        f"{accession}"
    )

    data = request_json(url)

    if not data:
        return None

    if isinstance(data, list):
        if not data:
            return None

        return data[0]

    return data


def download_alphafold_structure(
    gene_symbol: str,
    accession: str,
) -> Optional[Path]:

    output_file = (
        CACHE_DIR
        / f"{gene_symbol}_{accession}.cif"
    )

    if output_file.exists() and output_file.stat().st_size > 1000:
        return output_file

    metadata = get_alphafold_metadata(accession)

    if not metadata:
        logger.warning(
            "AlphaFold metadata unavailable for %s",
            gene_symbol,
        )
        return None

    cif_url = metadata.get("cifUrl")

    if not cif_url:
        logger.warning(
            "No CIF URL available for %s",
            gene_symbol,
        )
        return None

    try:
        response = SESSION.get(
            cif_url,
            timeout=REQUEST_TIMEOUT,
        )

        response.raise_for_status()

        if len(response.content) < 1000:
            logger.warning(
                "Downloaded CIF appears incomplete: %s",
                gene_symbol,
            )
            return None

        output_file.write_bytes(response.content)

        return output_file

    except requests.RequestException as exc:
        logger.warning(
            "AlphaFold download failed for %s: %s",
            gene_symbol,
            exc,
        )

    return None


# ============================================================
# 8. PARSE ALPHAFOLD STRUCTURE
# ============================================================

def extract_ca_coordinates(
    cif_path: Path,
) -> tuple[np.ndarray, np.ndarray]:

    parser = MMCIFParser(QUIET=True)

    structure = parser.get_structure(
        "protein",
        str(cif_path),
    )

    coordinates = []
    confidence = []

    for model in structure:

        for chain in model:

            for residue in chain:

                if not is_aa(residue, standard=False):
                    continue

                if "CA" not in residue:
                    continue

                atom = residue["CA"]

                coordinates.append(atom.coord)
                confidence.append(float(atom.bfactor))

            # Use the first chain to avoid combining
            # unrelated chains or biological assemblies.
            if coordinates:
                break

        if coordinates:
            break

    if not coordinates:
        raise ValueError(
            "No alpha-carbon coordinates found."
        )

    coordinates = np.asarray(
        coordinates,
        dtype=float,
    )

    confidence = np.asarray(
        confidence,
        dtype=float,
    )

    return coordinates, confidence


# ============================================================
# 9. 3D STRUCTURE VISUALIZATION
# ============================================================

def generate_3d_structure(
    gene: Gene,
    accession: str,
    cif_path: Path,
) -> bool:

    try:
        coordinates, confidence = extract_ca_coordinates(
            cif_path
        )

        if len(coordinates) > MAX_STRUCTURE_POINTS:

            indices = np.linspace(
                0,
                len(coordinates) - 1,
                MAX_STRUCTURE_POINTS,
                dtype=int,
            )

            coordinates = coordinates[indices]
            confidence = confidence[indices]

        x = coordinates[:, 0]
        y = coordinates[:, 1]
        z = coordinates[:, 2]

        fig = go.Figure()

        # Protein backbone.
        fig.add_trace(
            go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode="lines",
                line=dict(
                    color="#466A9F",
                    width=4,
                ),
                name="Protein backbone",
            )
        )

        # Confidence values represented by marker colour.
        fig.add_trace(
            go.Scatter3d(
                x=x,
                y=y,
                z=z,
                mode="markers",
                marker=dict(
                    size=2.5,
                    color=confidence,
                    colorscale="Viridis",
                    cmin=0,
                    cmax=100,
                    colorbar=dict(
                        title="pLDDT",
                    ),
                    opacity=0.9,
                ),
                name="Residue confidence",
            )
        )

        mean_confidence = float(
            np.mean(confidence)
        )

        fig.update_layout(
            title=dict(
                text=(
                    f"{gene.symbol} | AlphaFold predicted structure"
                    f"<br><sup>UniProt: {accession} | "
                    f"Mean pLDDT: {mean_confidence:.1f}</sup>"
                ),
                x=0.5,
            ),
            scene=dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False),
                aspectmode="data",
                bgcolor="white",
            ),
            paper_bgcolor="white",
            plot_bgcolor="white",
            margin=dict(
                l=0,
                r=0,
                b=0,
                t=80,
            ),
            showlegend=False,
        )

        output_file = (
            STRUCTURE_DIR
            / f"{gene.symbol}_alphafold_3d.png"
        )

        fig.write_image(
            str(output_file),
            width=1800,
            height=1400,
            scale=2,
        )

        return True

    except Exception as exc:
        logger.exception(
            "3D visualization failed for %s: %s",
            gene.symbol,
            exc,
        )

        return False


# ============================================================
# 10. CONCEPTUAL BIOLOGICAL ILLUSTRATIONS
# ============================================================

CATEGORY_COLORS = {
    "extracellular_matrix": "#4D86A8",
    "membrane_receptor": "#8C6BB1",
    "secreted_signalling": "#D18A45",
    "intracellular_signalling": "#4C9B78",
    "cytoskeleton": "#B75B66",
}


CATEGORY_LABELS = {
    "extracellular_matrix": "EXTRACELLULAR MATRIX",
    "membrane_receptor": "MEMBRANE RECEPTOR",
    "secreted_signalling": "SECRETED SIGNALLING",
    "intracellular_signalling": "INTRACELLULAR SIGNALLING",
    "cytoskeleton": "CYTOSKELETON",
}


def draw_cell_context(ax, category: str) -> None:
    """
    Draw a simplified conceptual cellular environment.

    This is not a molecularly accurate topology diagram.
    """

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.axis("off")

    # Extracellular region.
    ax.add_patch(
        FancyBboxPatch(
            (0.4, 3.7),
            9.2,
            2.7,
            boxstyle="round,pad=0.12",
            facecolor="#EEF3F7",
            edgecolor="#B8C7D5",
            linewidth=1.5,
        )
    )

    ax.text(
        0.7,
        6.05,
        "EXTRACELLULAR ENVIRONMENT",
        fontsize=8,
        color="#536B7C",
        weight="bold",
    )

    # Cell membrane.
    ax.plot(
        [0.5, 9.5],
        [3.45, 3.45],
        color="#9B6B9E",
        linewidth=5,
        alpha=0.8,
    )

    ax.plot(
        [0.5, 9.5],
        [3.2, 3.2],
        color="#CBAACB",
        linewidth=3,
        alpha=0.8,
    )

    # Intracellular region.
    ax.add_patch(
        FancyBboxPatch(
            (0.4, 0.35),
            9.2,
            2.65,
            boxstyle="round,pad=0.12",
            facecolor="#F4F2EB",
            edgecolor="#D5D0C1",
            linewidth=1.5,
        )
    )

    ax.text(
        0.7,
        0.6,
        "INTRACELLULAR SPACE",
        fontsize=8,
        color="#766E5E",
        weight="bold",
    )

    # Conceptual ECM fibres.
    for i in range(4):
        y = 4.35 + i * 0.35

        ax.plot(
            [1.0, 3.2, 5.0, 7.0, 9.0],
            [
                y,
                y + 0.25,
                y - 0.15,
                y + 0.2,
                y,
            ],
            color="#7CA4BB",
            linewidth=2,
            alpha=0.55,
        )

    # Conceptual intracellular network.
    for x1, y1, x2, y2 in [
        (2, 1.3, 4, 2.2),
        (4, 2.2, 6, 1.2),
        (6, 1.2, 8, 2.1),
        (2, 1.3, 6, 1.2),
    ]:

        ax.plot(
            [x1, x2],
            [y1, y2],
            color="#A7A08D",
            linewidth=1.5,
            alpha=0.65,
        )

    # Emphasize the relevant compartment.
    if category == "extracellular_matrix":

        ax.add_patch(
            Ellipse(
                (5, 5.1),
                2.3,
                1.0,
                facecolor="#4D86A8",
                alpha=0.22,
                edgecolor="#4D86A8",
                linewidth=2,
            )
        )

    elif category == "membrane_receptor":

        ax.add_patch(
            FancyBboxPatch(
                (4.4, 2.9),
                1.2,
                1.0,
                boxstyle="round,pad=0.1",
                facecolor="#8C6BB1",
                alpha=0.35,
                edgecolor="#8C6BB1",
                linewidth=2,
            )
        )

    elif category == "secreted_signalling":

        ax.add_patch(
            Circle(
                (5, 5.1),
                0.55,
                facecolor="#D18A45",
                alpha=0.35,
                edgecolor="#D18A45",
                linewidth=2,
            )
        )

    elif category == "intracellular_signalling":

        ax.add_patch(
            Circle(
                (5, 1.7),
                0.65,
                facecolor="#4C9B78",
                alpha=0.35,
                edgecolor="#4C9B78",
                linewidth=2,
            )
        )

    elif category == "cytoskeleton":

        ax.plot(
            [3.5, 5.0, 6.5],
            [1.1, 2.3, 1.1],
            color="#B75B66",
            linewidth=5,
            alpha=0.55,
        )


def generate_schematic(gene: Gene) -> bool:

    try:
        color = CATEGORY_COLORS[gene.category]

        fig, ax = plt.subplots(
            figsize=FIGURE_SIZE
        )

        fig.patch.set_facecolor("white")

        draw_cell_context(
            ax,
            gene.category,
        )

        # Protein representation.
        protein = Circle(
            (5, 4.95 if gene.category in [
                "extracellular_matrix",
                "secreted_signalling",
            ] else (
                3.45 if gene.category == "membrane_receptor"
                else 1.7
            )),
            0.48,
            facecolor=color,
            edgecolor="white",
            linewidth=2,
            alpha=0.95,
            zorder=10,
        )

        ax.add_patch(protein)

        # Protein label.
        protein_y = (
            4.95
            if gene.category in [
                "extracellular_matrix",
                "secreted_signalling",
            ]
            else (
                3.45
                if gene.category == "membrane_receptor"
                else 1.7
            )
        )

        ax.text(
            5,
            protein_y,
            gene.symbol,
            color="white",
            fontsize=9,
            weight="bold",
            ha="center",
            va="center",
            zorder=11,
        )

        # Biological role.
        ax.text(
            5,
            6.85,
            gene.symbol,
            ha="center",
            va="center",
            fontsize=23,
            weight="bold",
            color="#263746",
        )

        ax.text(
            5,
            6.48,
            CATEGORY_LABELS[gene.category],
            ha="center",
            va="center",
            fontsize=9,
            weight="bold",
            color=color,
        )

        ax.text(
            5,
            0.05,
            gene.role,
            ha="center",
            va="bottom",
            fontsize=10,
            color="#34495E",
            weight="bold",
        )

        fig.text(
            0.5,
            0.015,
            (
                "CONCEPTUAL BIOLOGICAL ILLUSTRATION | "
                "Not to scale | Not a molecular topology map"
            ),
            ha="center",
            fontsize=7,
            color="#777777",
        )

        output_file = (
            SCHEMATIC_DIR
            / f"{gene.symbol}_schematic.png"
        )

        save_figure(
            fig,
            output_file,
        )

        return True

    except Exception as exc:

        logger.exception(
            "Schematic generation failed for %s: %s",
            gene.symbol,
            exc,
        )

        return False


# ============================================================
# 11. AUDIT TABLE
# ============================================================

def write_audit(records: list[dict]) -> None:

    output_file = (
        AUDIT_DIR
        / "gene_protein_illustration_audit.csv"
    )

    if not records:
        logger.warning("No audit records available.")
        return

    with output_file.open(
        "w",
        newline="",
        encoding="utf-8-sig",
    ) as file:

        writer = csv.DictWriter(
            file,
            fieldnames=list(records[0].keys()),
        )

        writer.writeheader()
        writer.writerows(records)


# ============================================================
# 12. MAIN PIPELINE
# ============================================================

def main():

    logger.info("=" * 65)
    logger.info("TRANSCRIPTOME-X")
    logger.info("SCRIPT 44 - GENE AND PROTEIN ILLUSTRATIONS")
    logger.info("=" * 65)

    logger.info(
        "Genes to process: %d",
        len(GENES),
    )

    logger.info(
        "Output directory: %s",
        OUTPUT_ROOT,
    )

    records = []

    schematic_success = 0
    structure_success = 0

    for index, gene in enumerate(GENES, start=1):

        logger.info(
            "[%d/%d] Processing %s",
            index,
            len(GENES),
            gene.symbol,
        )

        record = {
            "gene": gene.symbol,
            "category": gene.category,
            "biological_role": gene.role,
            "schematic_status": "NOT_GENERATED",
            "uniprot_accession": "",
            "alphafold_status": "NOT_AVAILABLE",
            "structure_status": "NOT_GENERATED",
        }

        # ----------------------------------------------------
        # A. Schematic illustration
        # ----------------------------------------------------

        if generate_schematic(gene):

            record["schematic_status"] = "SUCCESS"
            schematic_success += 1

        else:

            record["schematic_status"] = "FAILED"

        # ----------------------------------------------------
        # B. UniProt accession
        # ----------------------------------------------------

        accession = get_uniprot_accession(
            gene.symbol
        )

        if accession:

            record["uniprot_accession"] = accession

            logger.info(
                "%s UniProt accession: %s",
                gene.symbol,
                accession,
            )

            # ------------------------------------------------
            # C. AlphaFold structure
            # ------------------------------------------------

            cif_path = download_alphafold_structure(
                gene.symbol,
                accession,
            )

            if cif_path:

                record["alphafold_status"] = "DOWNLOADED"

                if generate_3d_structure(
                    gene,
                    accession,
                    cif_path,
                ):

                    record["structure_status"] = "SUCCESS"
                    structure_success += 1

                else:

                    record["structure_status"] = "FAILED"

            else:

                record["alphafold_status"] = "UNAVAILABLE"

        else:

            record["alphafold_status"] = "UNIPROT_NOT_FOUND"

        records.append(record)

        write_audit(records)

        time.sleep(REQUEST_DELAY)

    # ========================================================
    # FINAL SUMMARY
    # ========================================================

    logger.info("=" * 65)
    logger.info("SCRIPT 44 - EXECUTION SUMMARY")
    logger.info("=" * 65)

    logger.info(
        "Total genes: %d",
        len(GENES),
    )

    logger.info(
        "Schematic illustrations generated: %d/%d",
        schematic_success,
        len(GENES),
    )

    logger.info(
        "AlphaFold 3D structures generated: %d/%d",
        structure_success,
        len(GENES),
    )

    logger.info(
        "Total images generated: %d",
        schematic_success + structure_success,
    )

    logger.info(
        "Schematic directory: %s",
        SCHEMATIC_DIR,
    )

    logger.info(
        "3D structure directory: %s",
        STRUCTURE_DIR,
    )

    logger.info(
        "Audit directory: %s",
        AUDIT_DIR,
    )

    logger.info("=" * 65)


if __name__ == "__main__":
    main()
