
"""
TRANSCRIPTOME-X
SCRIPT 47 - PYTHON PUBLICATION-QUALITY PROTEIN RENDERING

Generates three molecular representations for 23 proteins:

1. CARTOON:
   C-alpha backbone representation.

2. MOLECULAR SURFACE:
   Approximate van der Waals molecular surface calculated
   from heavy-atom coordinates.

3. SURFACE + BACKBONE:
   Molecular surface combined with the C-alpha trace.

Outputs:
- High-resolution PNG images
- Interactive Plotly HTML visualizations
- Structural rendering audit CSV
- Execution log

Input:
Script 46 structural audit CSV.

Original structures and Script 46 outputs are never modified.

Scientific note:
The backbone representation is a C-alpha trace, not a
DSSP-derived secondary-structure ribbon.
The molecular surface is an approximate van der Waals
surface, not an experimentally measured surface.
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

from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import plotly.graph_objects as go

from Bio.PDB import MMCIFParser, PDBParser

from scipy.spatial import cKDTree

from skimage.measure import marching_cubes


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

INPUT_AUDIT = (
    BASE_DIR
    / "final_3d_images"
    / "audit"
    / "protein_3d_image_audit.csv"
)

OUTPUT_DIR = (
    BASE_DIR
    / "python_publication_3d_rendering"
)

CARTOON_PNG = OUTPUT_DIR / "cartoon" / "png"
CARTOON_HTML = OUTPUT_DIR / "cartoon" / "html"

SURFACE_PNG = OUTPUT_DIR / "molecular_surface" / "png"
SURFACE_HTML = OUTPUT_DIR / "molecular_surface" / "html"

COMBINED_PNG = OUTPUT_DIR / "surface_cartoon" / "png"
COMBINED_HTML = OUTPUT_DIR / "surface_cartoon" / "html"

AUDIT_DIR = OUTPUT_DIR / "audit"

AUDIT_PATH = (
    AUDIT_DIR
    / "protein_rendering_audit.csv"
)

LOG_PATH = AUDIT_DIR / "script47.log"

# PNG quality.
PNG_DPI = 300
FIGURE_SIZE = (12, 10)

# Molecular surface.
GRID_SPACING = 0.85
SURFACE_PADDING = 5.0
MAX_GRID_DIMENSION = 180

# Rendering.
MAX_BACKBONE_POINTS = 12000
MAX_SURFACE_FACES = 100000

# HTML.
HTML_WIDTH = 1100
HTML_HEIGHT = 850

# Atom selection.
EXCLUDED_ELEMENTS = {
    "H",
    "D",
}

# Van der Waals radii in Angstrom.
VDW_RADII = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "P": 1.80,
    "F": 1.47,
    "CL": 1.75,
    "SE": 1.90,
    "MG": 1.73,
    "CA": 2.31,
    "ZN": 1.39,
    "FE": 1.56,
    "CU": 1.40,
    "MN": 1.73,
    "NA": 2.27,
    "K": 2.75,
}

EXPECTED_PROTEINS = 23


# ============================================================
# DIRECTORIES
# ============================================================

ALL_DIRECTORIES = [
    CARTOON_PNG,
    CARTOON_HTML,
    SURFACE_PNG,
    SURFACE_HTML,
    COMBINED_PNG,
    COMBINED_HTML,
    AUDIT_DIR,
]

for directory in ALL_DIRECTORIES:
    directory.mkdir(parents=True, exist_ok=True)


# ============================================================
# LOGGING
# ============================================================

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[
        logging.FileHandler(
            LOG_PATH,
            encoding="utf-8",
        ),
        logging.StreamHandler(sys.stdout),
    ],
)

logger = logging.getLogger("TRANSCRIPTOME-X-SCRIPT47")


# ============================================================
# INPUT AUDIT
# ============================================================

def load_input_records() -> list[dict]:

    if not INPUT_AUDIT.exists():
        raise FileNotFoundError(
            f"Script 46 audit not found: {INPUT_AUDIT}"
        )

    records = []

    with INPUT_AUDIT.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        reader = csv.DictReader(handle)

        required = {
            "gene",
            "uniprot_accession",
            "source",
            "structure_file",
            "status",
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                "Missing audit columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:

            gene = row["gene"].strip()
            source = row["source"].strip()
            structure_file = row["structure_file"].strip()

            if not gene or not structure_file:
                continue

            if row["status"].strip() not in {
                "SUCCESS",
                "PARTIAL",
            }:
                continue

            path = Path(structure_file)

            if not path.is_absolute():
                path = PROJECT_DIR / path

            if not path.is_file():
                logger.warning(
                    "%s: missing structure file: %s",
                    gene,
                    path,
                )
                continue

            records.append({
                "gene": gene,
                "accession": row["uniprot_accession"].strip(),
                "source": source,
                "structure_file": path,
            })

    genes = [record["gene"] for record in records]

    if len(genes) != len(set(genes)):
        raise ValueError(
            "Duplicate gene entries in Script 46 audit."
        )

    if not records:
        raise RuntimeError(
            "No valid structures found."
        )

    return records


# ============================================================
# STRUCTURE PARSING
# ============================================================

def parse_structure(path: Path, gene: str) -> dict:

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

    atom_coordinates = []
    atom_radii = []
    atom_elements = []

    backbone_chains = []
    backbone_bfactors = []

    residue_count = 0

    for model in structure:

        for chain in model:

            chain_coordinates = []
            chain_bfactors = []

            for residue in chain:

                if residue.id[0] != " ":
                    continue

                if not residue.has_id("CA"):
                    continue

                ca_atom = residue["CA"]

                chain_coordinates.append(
                    np.asarray(ca_atom.coord, dtype=float)
                )

                chain_bfactors.append(
                    float(ca_atom.get_bfactor())
                )

                residue_count += 1

                for atom in residue:

                    if atom.is_disordered():
                        atom = atom.selected_child

                    element = (
                        atom.element or ""
                    ).strip().upper()

                    if not element:
                        element = atom.get_name()[0].upper()

                    if element in EXCLUDED_ELEMENTS:
                        continue

                    coordinates = np.asarray(
                        atom.coord,
                        dtype=float,
                    )

                    if not np.all(np.isfinite(coordinates)):
                        continue

                    radius = VDW_RADII.get(
                        element,
                        1.70,
                    )

                    atom_coordinates.append(coordinates)
                    atom_radii.append(radius)
                    atom_elements.append(element)

            if len(chain_coordinates) >= 2:

                backbone_chains.append({
                    "chain": str(chain.id),
                    "coordinates": np.asarray(
                        chain_coordinates,
                        dtype=float,
                    ),
                    "bfactors": np.asarray(
                        chain_bfactors,
                        dtype=float,
                    ),
                })

    if not atom_coordinates:
        raise ValueError(
            "No valid heavy atoms found."
        )

    if not backbone_chains:
        raise ValueError(
            "No valid C-alpha backbone found."
        )

    return {
        "atoms": np.asarray(atom_coordinates, dtype=float),
        "radii": np.asarray(atom_radii, dtype=float),
        "elements": atom_elements,
        "chains": backbone_chains,
        "residues": residue_count,
    }


# ============================================================
# COORDINATE CENTERING
# ============================================================

def center_structure(data: dict) -> dict:

    center = np.mean(data["atoms"], axis=0)

    data["atoms"] = data["atoms"] - center

    for chain in data["chains"]:
        chain["coordinates"] = (
            chain["coordinates"] - center
        )

    return data


# ============================================================
# MOLECULAR SURFACE CALCULATION
# ============================================================

def calculate_surface(data: dict) -> tuple:

    atoms = data["atoms"]
    radii = data["radii"]

    minimum = np.min(
        atoms - radii[:, None],
        axis=0,
    ) - SURFACE_PADDING

    maximum = np.max(
        atoms + radii[:, None],
        axis=0,
    ) + SURFACE_PADDING

    spacing = GRID_SPACING

    dimensions = np.ceil(
        (maximum - minimum) / spacing
    ).astype(int) + 1

    largest_dimension = int(np.max(dimensions))

    if largest_dimension > MAX_GRID_DIMENSION:

        spacing *= (
            largest_dimension / MAX_GRID_DIMENSION
        )

        dimensions = np.ceil(
            (maximum - minimum) / spacing
        ).astype(int) + 1

    logger.info(
        "Surface grid: %s x %s x %s | spacing %.3f Å",
        dimensions[0],
        dimensions[1],
        dimensions[2],
        spacing,
    )

    grid = np.full(
        tuple(dimensions),
        20.0,
        dtype=np.float32,
    )

    tree = cKDTree(atoms)

    # Calculate the minimum distance to atom surfaces.
    # Process in chunks to control memory consumption.

    nx, ny, nz = dimensions

    chunk_size = 12

    x_values = (
        minimum[0]
        + np.arange(nx) * spacing
    )

    y_values = (
        minimum[1]
        + np.arange(ny) * spacing
    )

    z_values = (
        minimum[2]
        + np.arange(nz) * spacing
    )

    for x_start in range(0, nx, chunk_size):

        x_chunk = x_values[
            x_start:x_start + chunk_size
        ]

        xx, yy, zz = np.meshgrid(
            x_chunk,
            y_values,
            z_values,
            indexing="ij",
        )

        points = np.column_stack([
            xx.ravel(),
            yy.ravel(),
            zz.ravel(),
        ])

        distances, indices = tree.query(
            points,
            k=1,
            workers=-1,
        )

        # Signed distance to the nearest atom's
        # van der Waals sphere.
        values = (
            distances - radii[indices]
        )

        grid[
            x_start:x_start + len(x_chunk),
            :,
            :,
        ] = values.reshape(
            len(x_chunk),
            ny,
            nz,
        ).astype(np.float32)

    vertices, faces, normals, values = marching_cubes(
        grid,
        level=0.0,
        spacing=(spacing, spacing, spacing),
    )

    vertices += minimum

    if len(faces) == 0:
        raise ValueError(
            "Molecular surface generation returned no faces."
        )

    return vertices, faces, normals


# ============================================================
# SURFACE MESH DOWNSAMPLING
# ============================================================

def sample_faces(
    vertices: np.ndarray,
    faces: np.ndarray,
    maximum: int,
) -> np.ndarray:

    if len(faces) <= maximum:
        return faces

    indices = np.linspace(
        0,
        len(faces) - 1,
        maximum,
        dtype=int,
    )

    return faces[indices]


# ============================================================
# BACKBONE PREPARATION
# ============================================================

def sample_backbone(
    chains: list[dict],
) -> list[dict]:

    result = []

    total_points = sum(
        len(chain["coordinates"])
        for chain in chains
    )

    maximum_per_chain = max(
        2,
        MAX_BACKBONE_POINTS // max(1, len(chains)),
    )

    for chain in chains:

        coordinates = chain["coordinates"]
        bfactors = chain["bfactors"]

        if len(coordinates) > maximum_per_chain:

            indices = np.linspace(
                0,
                len(coordinates) - 1,
                maximum_per_chain,
                dtype=int,
            )

            coordinates = coordinates[indices]
            bfactors = bfactors[indices]

        result.append({
            "chain": chain["chain"],
            "coordinates": coordinates,
            "bfactors": bfactors,
        })

    return result


# ============================================================
# MATPLOTLIB HELPERS
# ============================================================

def set_3d_axes(ax, vertices: np.ndarray) -> None:

    minimum = np.min(vertices, axis=0)
    maximum = np.max(vertices, axis=0)

    center = (minimum + maximum) / 2.0

    extent = float(
        np.max(maximum - minimum)
    )

    half = max(extent / 2.0, 1.0)

    ax.set_xlim(center[0] - half, center[0] + half)
    ax.set_ylim(center[1] - half, center[1] + half)
    ax.set_zlim(center[2] - half, center[2] + half)

    ax.set_box_aspect((1, 1, 1))

    ax.set_axis_off()

    ax.view_init(
        elev=24,
        azim=42,
    )


def add_surface_matplotlib(
    ax,
    vertices: np.ndarray,
    faces: np.ndarray,
    color: str = "lightsteelblue",
    alpha: float = 1.0,
) -> None:

    selected_faces = sample_faces(
        vertices,
        faces,
        MAX_SURFACE_FACES,
    )

    triangles = vertices[selected_faces]

    collection = Poly3DCollection(
        triangles,
        facecolor=color,
        edgecolor="none",
        alpha=alpha,
        linewidth=0,
    )

    collection.set_zsort("average")

    ax.add_collection3d(collection)


def add_backbone_matplotlib(
    ax,
    chains: list[dict],
) -> None:

    for index, chain in enumerate(chains):

        coordinates = chain["coordinates"]

        ax.plot(
            coordinates[:, 0],
            coordinates[:, 1],
            coordinates[:, 2],
            linewidth=1.5,
            alpha=0.95,
            color=(
                "firebrick"
                if index % 3 == 0
                else "royalblue"
                if index % 3 == 1
                else "darkorange"
            ),
        )


def save_matplotlib_figure(
    fig,
    output_path: Path,
) -> None:

    fig.savefig(
        output_path,
        dpi=PNG_DPI,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(fig)


# ============================================================
# CARTOON PNG
# ============================================================

def generate_cartoon_png(
    gene: str,
    chains: list[dict],
    output_path: Path,
) -> None:

    fig = plt.figure(
        figsize=FIGURE_SIZE,
    )

    ax = fig.add_subplot(
        111,
        projection="3d",
    )

    all_coordinates = np.vstack([
        chain["coordinates"]
        for chain in chains
    ])

    add_backbone_matplotlib(
        ax,
        chains,
    )

    set_3d_axes(
        ax,
        all_coordinates,
    )

    ax.set_title(
        f"{gene} | C-alpha backbone",
        fontsize=17,
        pad=18,
    )

    save_matplotlib_figure(
        fig,
        output_path,
    )


# ============================================================
# SURFACE PNG
# ============================================================

def generate_surface_png(
    gene: str,
    vertices: np.ndarray,
    faces: np.ndarray,
    output_path: Path,
) -> None:

    fig = plt.figure(
        figsize=FIGURE_SIZE,
    )

    ax = fig.add_subplot(
        111,
        projection="3d",
    )

    add_surface_matplotlib(
        ax,
        vertices,
        faces,
        color="lightsteelblue",
        alpha=1.0,
    )

    set_3d_axes(
        ax,
        vertices,
    )

    ax.set_title(
        f"{gene} | Molecular surface",
        fontsize=17,
        pad=18,
    )

    save_matplotlib_figure(
        fig,
        output_path,
    )


# ============================================================
# COMBINED PNG
# ============================================================

def generate_combined_png(
    gene: str,
    vertices: np.ndarray,
    faces: np.ndarray,
    chains: list[dict],
    output_path: Path,
) -> None:

    fig = plt.figure(
        figsize=FIGURE_SIZE,
    )

    ax = fig.add_subplot(
        111,
        projection="3d",
    )

    add_surface_matplotlib(
        ax,
        vertices,
        faces,
        color="lightsteelblue",
        alpha=0.32,
    )

    add_backbone_matplotlib(
        ax,
        chains,
    )

    set_3d_axes(
        ax,
        vertices,
    )

    ax.set_title(
        f"{gene} | Surface + C-alpha backbone",
        fontsize=17,
        pad=18,
    )

    save_matplotlib_figure(
        fig,
        output_path,
    )


# ============================================================
# PLOTLY SURFACE TRACE
# ============================================================

def make_surface_trace(
    vertices: np.ndarray,
    faces: np.ndarray,
    opacity: float = 1.0,
    color: str = "lightsteelblue",
) -> go.Mesh3d:

    selected_faces = sample_faces(
        vertices,
        faces,
        MAX_SURFACE_FACES,
    )

    return go.Mesh3d(
        x=vertices[:, 0],
        y=vertices[:, 1],
        z=vertices[:, 2],
        i=selected_faces[:, 0],
        j=selected_faces[:, 1],
        k=selected_faces[:, 2],
        color=color,
        opacity=opacity,
        flatshading=False,
        lighting={
            "ambient": 0.45,
            "diffuse": 0.75,
            "specular": 0.35,
            "roughness": 0.55,
            "fresnel": 0.15,
        },
        name="Molecular surface",
        hoverinfo="skip",
    )


# ============================================================
# PLOTLY BACKBONE TRACE
# ============================================================

def make_backbone_traces(
    chains: list[dict],
) -> list[go.Scatter3d]:

    traces = []

    palette = [
        "#B64040",
        "#2867A5",
        "#D58A2B",
        "#36876D",
    ]

    for index, chain in enumerate(chains):

        coordinates = chain["coordinates"]

        traces.append(
            go.Scatter3d(
                x=coordinates[:, 0],
                y=coordinates[:, 1],
                z=coordinates[:, 2],
                mode="lines",
                name=f"Chain {chain['chain']}",
                line={
                    "color": palette[index % len(palette)],
                    "width": 5,
                },
                hoverinfo="skip",
            )
        )

    return traces


# ============================================================
# PLOTLY LAYOUT
# ============================================================

def configure_plotly(
    fig: go.Figure,
    title: str,
) -> None:

    fig.update_layout(
        title={
            "text": title,
            "x": 0.5,
            "xanchor": "center",
        },
        width=HTML_WIDTH,
        height=HTML_HEIGHT,
        template="plotly_white",
        paper_bgcolor="white",
        plot_bgcolor="white",
        margin={
            "l": 0,
            "r": 0,
            "t": 65,
            "b": 0,
        },
        scene={
            "aspectmode": "data",
            "xaxis": {"visible": False},
            "yaxis": {"visible": False},
            "zaxis": {"visible": False},
            "bgcolor": "white",
        },
        legend={
            "orientation": "h",
            "y": 0.02,
        },
    )


def save_plotly(
    fig: go.Figure,
    output_path: Path,
    title: str,
) -> None:

    configure_plotly(
        fig,
        title,
    )

    fig.write_html(
        str(output_path),
        include_plotlyjs=True,
        full_html=True,
    )


# ============================================================
# CARTOON HTML
# ============================================================

def generate_cartoon_html(
    gene: str,
    chains: list[dict],
    output_path: Path,
) -> None:

    fig = go.Figure()

    for trace in make_backbone_traces(chains):
        fig.add_trace(trace)

    save_plotly(
        fig,
        output_path,
        f"{gene} | C-alpha backbone",
    )


# ============================================================
# SURFACE HTML
# ============================================================

def generate_surface_html(
    gene: str,
    vertices: np.ndarray,
    faces: np.ndarray,
    output_path: Path,
) -> None:

    fig = go.Figure()

    fig.add_trace(
        make_surface_trace(
            vertices,
            faces,
        )
    )

    save_plotly(
        fig,
        output_path,
        f"{gene} | Molecular surface",
    )


# ============================================================
# COMBINED HTML
# ============================================================

def generate_combined_html(
    gene: str,
    vertices: np.ndarray,
    faces: np.ndarray,
    chains: list[dict],
    output_path: Path,
) -> None:

    fig = go.Figure()

    fig.add_trace(
        make_surface_trace(
            vertices,
            faces,
            opacity=0.38,
        )
    )

    for trace in make_backbone_traces(chains):
        fig.add_trace(trace)

    save_plotly(
        fig,
        output_path,
        f"{gene} | Surface + C-alpha backbone",
    )


# ============================================================
# OUTPUT VALIDATION
# ============================================================

def valid_file(path: Path) -> bool:

    return (
        path.is_file()
        and path.stat().st_size > 1000
    )


# ============================================================
# SINGLE PROTEIN PIPELINE
# ============================================================

def process_protein(record: dict) -> dict:

    gene = record["gene"]

    result = {
        "gene": gene,
        "uniprot_accession": record["accession"],
        "source": record["source"],
        "structure_file": str(record["structure_file"]),
        "heavy_atoms": 0,
        "ca_residues": 0,
        "surface_vertices": 0,
        "surface_faces": 0,
        "cartoon_png": "NOT_GENERATED",
        "cartoon_html": "NOT_GENERATED",
        "surface_png": "NOT_GENERATED",
        "surface_html": "NOT_GENERATED",
        "combined_png": "NOT_GENERATED",
        "combined_html": "NOT_GENERATED",
        "status": "PENDING",
        "error": "",
    }

    logger.info(
        "Processing %s",
        gene,
    )

    try:

        data = parse_structure(
            record["structure_file"],
            gene,
        )

        data = center_structure(data)

        chains = sample_backbone(
            data["chains"]
        )

        result["heavy_atoms"] = len(data["atoms"])

        result["ca_residues"] = data["residues"]

        logger.info(
            "%s: %s heavy atoms | %s residues",
            gene,
            result["heavy_atoms"],
            result["ca_residues"],
        )

        # Calculate molecular surface once.
        vertices, faces, normals = calculate_surface(data)

        result["surface_vertices"] = len(vertices)
        result["surface_faces"] = len(faces)

        # Output paths.
        paths = {
            "cartoon_png": CARTOON_PNG / f"{gene}_cartoon.png",
            "cartoon_html": CARTOON_HTML / f"{gene}_cartoon.html",
            "surface_png": SURFACE_PNG / f"{gene}_surface.png",
            "surface_html": SURFACE_HTML / f"{gene}_surface.html",
            "combined_png": COMBINED_PNG / f"{gene}_surface_cartoon.png",
            "combined_html": COMBINED_HTML / f"{gene}_surface_cartoon.html",
        }

        # 1. Cartoon.
        generate_cartoon_png(
            gene,
            chains,
            paths["cartoon_png"],
        )

        generate_cartoon_html(
            gene,
            chains,
            paths["cartoon_html"],
        )

        # 2. Molecular surface.
        generate_surface_png(
            gene,
            vertices,
            faces,
            paths["surface_png"],
        )

        generate_surface_html(
            gene,
            vertices,
            faces,
            paths["surface_html"],
        )

        # 3. Combined.
        generate_combined_png(
            gene,
            vertices,
            faces,
            chains,
            paths["combined_png"],
        )

        generate_combined_html(
            gene,
            vertices,
            faces,
            chains,
            paths["combined_html"],
        )

        for field, path in paths.items():

            result[field] = (
                "SUCCESS"
                if valid_file(path)
                else "FAILED"
            )

        statuses = [
            result[field]
            for field in paths
        ]

        if all(status == "SUCCESS" for status in statuses):
            result["status"] = "SUCCESS"

        elif any(status == "SUCCESS" for status in statuses):
            result["status"] = "PARTIAL"

        else:
            result["status"] = "FAILED"

        logger.info(
            "%s completed: %s",
            gene,
            result["status"],
        )

    except Exception as exc:

        result["status"] = "ERROR"
        result["error"] = str(exc)

        logger.exception(
            "%s rendering failed.",
            gene,
        )

    return result


# ============================================================
# SAVE AUDIT
# ============================================================

def save_audit(records: list[dict]) -> None:

    if not records:
        return

    with AUDIT_PATH.open(
        "w",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=list(records[0].keys()),
        )

        writer.writeheader()
        writer.writerows(records)


# ============================================================
# MAIN
# ============================================================

def main():

    logger.info("=" * 75)
    logger.info("TRANSCRIPTOME-X - SCRIPT 47")
    logger.info("PYTHON PUBLICATION-QUALITY PROTEIN RENDERING")
    logger.info("=" * 75)

    logger.info(
        "Reading Script 46 structural audit..."
    )

    records = load_input_records()

    logger.info(
        "Structures available: %s/%s",
        len(records),
        EXPECTED_PROTEINS,
    )

    audit_records = []

    for index, record in enumerate(
        records,
        start=1,
    ):

        logger.info(
            "\n[%s/%s] %s",
            index,
            len(records),
            record["gene"],
        )

        result = process_protein(record)

        audit_records.append(result)

        # Save progressively to preserve completed results
        # if a later protein fails.
        save_audit(audit_records)

    # ========================================================
    # FINAL SUMMARY
    # ========================================================

    total = len(audit_records)

    def count_success(field: str) -> int:

        return sum(
            record[field] == "SUCCESS"
            for record in audit_records
        )

    complete = sum(
        record["status"] == "SUCCESS"
        for record in audit_records
    )

    print("\n" + "=" * 75)
    print("SCRIPT 47 - FINAL SUMMARY")
    print("=" * 75)

    print(f"Proteins processed: {total}")

    print(
        f"Cartoon PNG: "
        f"{count_success('cartoon_png')}/{total}"
    )

    print(
        f"Cartoon HTML: "
        f"{count_success('cartoon_html')}/{total}"
    )

    print(
        f"Molecular Surface PNG: "
        f"{count_success('surface_png')}/{total}"
    )

    print(
        f"Molecular Surface HTML: "
        f"{count_success('surface_html')}/{total}"
    )

    print(
        f"Combined PNG: "
        f"{count_success('combined_png')}/{total}"
    )

    print(
        f"Combined HTML: "
        f"{count_success('combined_html')}/{total}"
    )

    print(
        f"\nCompletely successful proteins: "
        f"{complete}/{total}"
    )

    print("\nRESULTS BY PROTEIN")

    for record in audit_records:

        print(
            f"{record['gene']:10s} | "
            f"{record['status']:10s} | "
            f"Atoms: {record['heavy_atoms']:7d} | "
            f"Residues: {record['ca_residues']:5d}"
        )

    print("\nOUTPUT DIRECTORIES")

    print(f"Cartoon PNG: {CARTOON_PNG}")
    print(f"Cartoon HTML: {CARTOON_HTML}")

    print(f"Surface PNG: {SURFACE_PNG}")
    print(f"Surface HTML: {SURFACE_HTML}")

    print(f"Combined PNG: {COMBINED_PNG}")
    print(f"Combined HTML: {COMBINED_HTML}")

    print(f"Audit: {AUDIT_PATH}")

    print("=" * 75)


if __name__ == "__main__":
    main()
