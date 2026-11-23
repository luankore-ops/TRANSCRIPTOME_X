
# ============================================================
# TRANSCRIPTOME-X
# SCRIPT 48 - PyMOL PROTEIN STRUCTURAL RENDERING
# ============================================================
#
# Generates publication-quality molecular visualizations:
#   1. Cartoon / Ribbon
#   2. Molecular Surface
#   3. Surface + Cartoon
#
# Input:
#   results/gene_protein_illustrations/final_3d_images/
#       audit/protein_3d_image_audit.csv
#
# Output:
#   results/gene_protein_illustrations/pymol_3d_rendering/
#
# Requirements:
#   pymol-open-source
#   biopython
#   pandas
#   dssp (installed through Conda)
#
# Run:
#   python 48_pymol_protein_rendering.py
#
# ============================================================

from __future__ import annotations

import csv
import logging
import os
import re
import sys
import traceback

from datetime import datetime
from pathlib import Path

import pandas as pd

try:
    import pymol2
except ImportError as exc:
    raise SystemExit(
        "ERROR: PyMOL is not available in this Python environment.\n"
        "Activate the transcriptome_pymol Conda environment first."
    ) from exc


# ============================================================
# 1. PROJECT CONFIGURATION
# ============================================================

PROJECT_ROOT = Path(
    r"C:\Users\luank\Dropbox\PC\Desktop\PYTHON\TRANSCRIPTOME_X"
)

INPUT_AUDIT = (
    PROJECT_ROOT
    / "results"
    / "gene_protein_illustrations"
    / "final_3d_images"
    / "audit"
    / "protein_3d_image_audit.csv"
)

OUTPUT_ROOT = (
    PROJECT_ROOT
    / "results"
    / "gene_protein_illustrations"
    / "pymol_3d_rendering"
)

OUTPUT_CARTOON = OUTPUT_ROOT / "cartoon" / "png"
OUTPUT_SURFACE = OUTPUT_ROOT / "molecular_surface" / "png"
OUTPUT_COMBINED = OUTPUT_ROOT / "surface_cartoon" / "png"

OUTPUT_SESSIONS = OUTPUT_ROOT / "sessions"
OUTPUT_SCRIPTS = OUTPUT_ROOT / "pymol_scripts"
OUTPUT_AUDIT = OUTPUT_ROOT / "audit"

IMAGE_WIDTH = 2400
IMAGE_HEIGHT = 1800

RAY_QUALITY = 2
ANTIALIAS = 2

# Surface transparency in the combined representation.
COMBINED_SURFACE_TRANSPARENCY = 0.28

# Use a white background for consistent scientific figures.
BACKGROUND_COLOR = "white"

# Protein secondary structure colors.
HELIX_COLOR = "salmon"
SHEET_COLOR = "marine"
LOOP_COLOR = "gray80"

# Default protein color.
DEFAULT_COLOR = "gray70"

# ============================================================
# 2. LOGGING
# ============================================================

OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

for folder in (
    OUTPUT_CARTOON,
    OUTPUT_SURFACE,
    OUTPUT_COMBINED,
    OUTPUT_SESSIONS,
    OUTPUT_SCRIPTS,
    OUTPUT_AUDIT,
):
    folder.mkdir(parents=True, exist_ok=True)

LOG_FILE = OUTPUT_AUDIT / "script48.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE, encoding="utf-8"),
        logging.StreamHandler(sys.stdout),
    ],
)

logger = logging.getLogger("TRANSCRIPTOME_X_SCRIPT48")


# ============================================================
# 3. HELPERS
# ============================================================

def safe_filename(value: str) -> str:
    """Convert a protein identifier into a safe filename."""

    value = str(value).strip()

    value = re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        value,
    )

    return value.strip("_") or "UNKNOWN"


def resolve_structure_path(value: str) -> Path:
    """Resolve structure paths from the Script 46 audit."""

    raw_path = str(value).strip().strip('"')

    if not raw_path or raw_path.lower() in ("nan", "none"):
        raise ValueError("Missing structure_file path.")

    path = Path(raw_path)

    if not path.is_absolute():
        path = PROJECT_ROOT / path

    path = path.resolve()

    if not path.exists():
        raise FileNotFoundError(
            f"Structure file not found: {path}"
        )

    if path.suffix.lower() not in (".pdb", ".ent", ".cif", ".mmcif"):
        raise ValueError(
            f"Unsupported structure format: {path.suffix}"
        )

    return path


def configure_scene(cmd):
    """Apply consistent rendering settings."""

    cmd.bg_color(BACKGROUND_COLOR)

    cmd.set("orthoscopic", 1)
    cmd.set("depth_cue", 0)
    cmd.set("ray_opaque_background", 1)

    cmd.set("antialias", ANTIALIAS)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ray_trace_gain", 0.08)
    cmd.set("ray_shadows", 1)

    cmd.set("specular", 0.25)
    cmd.set("shininess", 30)

    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("cartoon_flat_sheets", 1)
    cmd.set("cartoon_sampling", 14)
    cmd.set("cartoon_transparency", 0)

    cmd.set("surface_quality", 1)
    cmd.set("surface_smooth_edges", 1)

    cmd.set("transparency_mode", 2)

    cmd.set("two_sided_lighting", 1)
    cmd.set("light_count", 2)

    cmd.set("ray_trace_fog", 0)
    cmd.set("ray_trace_disco_factor", 1)

    cmd.set("label_size", 18)


def prepare_protein(cmd, structure_path: Path):
    """
    Load the structure and prepare its protein representation.

    PyMOL's DSS assignment is used to identify helices,
    sheets and loops for the cartoon representation.
    """

    cmd.reinitialize()

    configure_scene(cmd)

    cmd.load(
        str(structure_path),
        "protein",
    )

    cmd.remove("solvent")

    # Keep protein atoms only.
    cmd.remove("not polymer.protein")

    cmd.hide("everything", "all")

    # Assign secondary structure.
    cmd.dss("protein")

    # Base color.
    cmd.color(
        DEFAULT_COLOR,
        "protein",
    )

    # Secondary structure coloring.
    cmd.color(
        HELIX_COLOR,
        "protein and ss H",
    )

    cmd.color(
        SHEET_COLOR,
        "protein and ss S",
    )

    cmd.color(
        LOOP_COLOR,
        "protein and ss L",
    )

    # Improve molecular appearance.
    cmd.orient("protein")
    cmd.zoom("protein", 1.12)

    return True


def save_png(cmd, output_path: Path):
    """Render and save a high-resolution PNG."""

    cmd.png(
        str(output_path),
        width=IMAGE_WIDTH,
        height=IMAGE_HEIGHT,
        dpi=300,
        ray=1,
    )

    if not output_path.exists():
        raise RuntimeError(
            f"PNG generation failed: {output_path}"
        )

    if output_path.stat().st_size == 0:
        raise RuntimeError(
            f"Generated PNG is empty: {output_path}"
        )


# ============================================================
# 4. REPRESENTATION 1 - CARTOON / RIBBON
# ============================================================

def render_cartoon(cmd, output_path: Path):
    """Render a secondary-structure-aware cartoon."""

    cmd.hide("everything", "all")

    cmd.show(
        "cartoon",
        "protein",
    )

    cmd.set(
        "cartoon_transparency",
        0,
    )

    cmd.set(
        "cartoon_flat_sheets",
        1,
    )

    cmd.set(
        "cartoon_fancy_helices",
        1,
    )

    cmd.set(
        "cartoon_smooth_loops",
        1,
    )

    cmd.orient("protein")
    cmd.zoom("protein", 1.12)

    save_png(cmd, output_path)


# ============================================================
# 5. REPRESENTATION 2 - MOLECULAR SURFACE
# ============================================================

def render_surface(cmd, output_path: Path):
    """Render the molecular surface."""

    cmd.hide("everything", "all")

    cmd.show(
        "surface",
        "protein",
    )

    cmd.set(
        "transparency",
        0,
    )

    cmd.set(
        "surface_quality",
        1,
    )

    # Preserve secondary-structure color assignments.
    cmd.orient("protein")
    cmd.zoom("protein", 1.12)

    save_png(cmd, output_path)


# ============================================================
# 6. REPRESENTATION 3 - SURFACE + CARTOON
# ============================================================

def render_combined(cmd, output_path: Path):
    """Render a partially transparent surface with cartoon."""

    cmd.hide("everything", "all")

    # Molecular surface.
    cmd.show(
        "surface",
        "protein",
    )

    cmd.set(
        "transparency",
        COMBINED_SURFACE_TRANSPARENCY,
    )

    # Internal cartoon.
    cmd.show(
        "cartoon",
        "protein",
    )

    cmd.set(
        "cartoon_transparency",
        0,
    )

    cmd.set(
        "cartoon_fancy_helices",
        1,
    )

    cmd.set(
        "cartoon_smooth_loops",
        1,
    )

    cmd.orient("protein")
    cmd.zoom("protein", 1.12)

    save_png(cmd, output_path)


# ============================================================
# 7. GENERATE REPRODUCIBLE PYMOL SCRIPT
# ============================================================

def write_pymol_script(
    gene: str,
    structure_path: Path,
    script_path: Path,
):
    """Write a standalone PyMOL script for the combined view."""

    structure = str(structure_path).replace("\\", "/")

    commands = [
        "# TRANSCRIPTOME-X - PyMOL rendering",
        f"# Protein: {gene}",
        "",
        f'load "{structure}", protein',
        "remove solvent",
        "remove not polymer.protein",
        "hide everything, all",
        "dss protein",
        "",
        f"color {DEFAULT_COLOR}, protein",
        f"color {HELIX_COLOR}, protein and ss H",
        f"color {SHEET_COLOR}, protein and ss S",
        f"color {LOOP_COLOR}, protein and ss L",
        "",
        "bg_color white",
        "set orthoscopic, 1",
        "set depth_cue, 0",
        "set antialias, 2",
        "set ray_opaque_background, 1",
        "",
        "show surface, protein",
        f"set transparency, {COMBINED_SURFACE_TRANSPARENCY}",
        "show cartoon, protein",
        "set cartoon_fancy_helices, 1",
        "set cartoon_smooth_loops, 1",
        "set cartoon_flat_sheets, 1",
        "",
        "orient protein",
        "zoom protein, 1.12",
        "",
        "# Render manually using:",
        "# ray 2400, 1800",
        "# png output.png, dpi=300",
    ]

    script_path.write_text(
        "\n".join(commands),
        encoding="utf-8",
    )


# ============================================================
# 8. PROCESS ONE PROTEIN
# ============================================================

def process_protein(
    pymol_instance,
    record: dict,
) -> dict:

    gene = str(record["gene"]).strip()

    accession = str(
        record.get("uniprot_accession", "")
    ).strip()

    source = str(
        record.get("source", "")
    ).strip()

    structure_path = resolve_structure_path(
        record["structure_file"]
    )

    identifier = safe_filename(gene)

    logger.info("=" * 65)
    logger.info("PROCESSING: %s", gene)
    logger.info("SOURCE: %s", source)
    logger.info("STRUCTURE: %s", structure_path)

    result = {
        "gene": gene,
        "uniprot_accession": accession,
        "source": source,
        "structure_file": str(structure_path),
        "cartoon_status": "NOT_GENERATED",
        "surface_status": "NOT_GENERATED",
        "combined_status": "NOT_GENERATED",
        "session_status": "NOT_GENERATED",
        "script_status": "NOT_GENERATED",
        "status": "FAILED",
        "error": "",
    }

    cmd = pymol_instance.cmd

    try:
        prepare_protein(
            cmd,
            structure_path,
        )

        # ----------------------------------------------------
        # CARTOON
        # ----------------------------------------------------

        cartoon_path = (
            OUTPUT_CARTOON
            / f"{identifier}_cartoon.png"
        )

        render_cartoon(
            cmd,
            cartoon_path,
        )

        result["cartoon_status"] = "SUCCESS"

        logger.info(
            "Cartoon generated: %s",
            cartoon_path.name,
        )

        # ----------------------------------------------------
        # MOLECULAR SURFACE
        # ----------------------------------------------------

        surface_path = (
            OUTPUT_SURFACE
            / f"{identifier}_molecular_surface.png"
        )

        render_surface(
            cmd,
            surface_path,
        )

        result["surface_status"] = "SUCCESS"

        logger.info(
            "Surface generated: %s",
            surface_path.name,
        )

        # ----------------------------------------------------
        # SURFACE + CARTOON
        # ----------------------------------------------------

        combined_path = (
            OUTPUT_COMBINED
            / f"{identifier}_surface_cartoon.png"
        )

        render_combined(
            cmd,
            combined_path,
        )

        result["combined_status"] = "SUCCESS"

        logger.info(
            "Combined rendering generated: %s",
            combined_path.name,
        )

        # ----------------------------------------------------
        # PYMOL SESSION
        # ----------------------------------------------------

        session_path = (
            OUTPUT_SESSIONS
            / f"{identifier}.pse"
        )

        cmd.save(
            str(session_path)
        )

        if session_path.exists():
            result["session_status"] = "SUCCESS"

        # ----------------------------------------------------
        # REPRODUCIBLE PYMOL SCRIPT
        # ----------------------------------------------------

        script_path = (
            OUTPUT_SCRIPTS
            / f"{identifier}.pml"
        )

        write_pymol_script(
            gene,
            structure_path,
            script_path,
        )

        result["script_status"] = "SUCCESS"

        # ----------------------------------------------------
        # FINAL STATUS
        # ----------------------------------------------------

        statuses = [
            result["cartoon_status"],
            result["surface_status"],
            result["combined_status"],
            result["session_status"],
            result["script_status"],
        ]

        result["status"] = (
            "SUCCESS"
            if all(x == "SUCCESS" for x in statuses)
            else "PARTIAL"
        )

    except Exception as exc:

        result["error"] = str(exc)

        logger.error(
            "FAILED: %s | %s",
            gene,
            exc,
        )

        logger.debug(
            traceback.format_exc()
        )

    return result


# ============================================================
# 9. MAIN
# ============================================================

def main():

    logger.info("=" * 65)
    logger.info("TRANSCRIPTOME-X")
    logger.info("SCRIPT 48 - PYMOL PROTEIN STRUCTURAL RENDERING")
    logger.info("=" * 65)

    logger.info("Project: %s", PROJECT_ROOT)
    logger.info("Input audit: %s", INPUT_AUDIT)
    logger.info("Output: %s", OUTPUT_ROOT)

    if not PROJECT_ROOT.exists():
        raise FileNotFoundError(
            f"Project directory not found: {PROJECT_ROOT}"
        )

    if not INPUT_AUDIT.exists():
        raise FileNotFoundError(
            f"Input audit not found: {INPUT_AUDIT}"
        )

    # --------------------------------------------------------
    # LOAD INPUT AUDIT
    # --------------------------------------------------------

    df = pd.read_csv(
        INPUT_AUDIT,
        dtype=str,
    )

    required_columns = [
        "gene",
        "uniprot_accession",
        "source",
        "structure_file",
        "status",
    ]

    missing = [
        column
        for column in required_columns
        if column not in df.columns
    ]

    if missing:
        raise ValueError(
            f"Missing required audit columns: {missing}"
        )

    df = df[
        df["status"].str.upper() == "SUCCESS"
    ].copy()

    df = df.drop_duplicates(
        subset=["gene"],
        keep="first",
    )

    logger.info(
        "Proteins available for rendering: %d",
        len(df),
    )

    if df.empty:
        raise RuntimeError(
            "No proteins with SUCCESS status were found."
        )

    # --------------------------------------------------------
    # INITIALIZE PYMOL
    # --------------------------------------------------------

    pymol_instance = pymol2.PyMOL()

    results = []

    try:

        pymol_instance.start()

        logger.info("PyMOL initialized successfully.")

        # ----------------------------------------------------
        # PROCESS PROTEINS
        # ----------------------------------------------------

        for index, (_, row) in enumerate(
            df.iterrows(),
            start=1,
        ):

            logger.info(
                "Protein %d/%d",
                index,
                len(df),
            )

            record = row.to_dict()

            result = process_protein(
                pymol_instance,
                record,
            )

            results.append(result)

    finally:

        try:
            pymol_instance.stop()
        except Exception:
            logger.exception(
                "Error while shutting down PyMOL."
            )

    # --------------------------------------------------------
    # SAVE AUDIT
    # --------------------------------------------------------

    audit_df = pd.DataFrame(results)

    audit_path = (
        OUTPUT_AUDIT
        / "protein_pymol_rendering_audit.csv"
    )

    audit_df.to_csv(
        audit_path,
        index=False,
        encoding="utf-8-sig",
    )

    # --------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------

    total = len(results)

    cartoon_success = sum(
        r["cartoon_status"] == "SUCCESS"
        for r in results
    )

    surface_success = sum(
        r["surface_status"] == "SUCCESS"
        for r in results
    )

    combined_success = sum(
        r["combined_status"] == "SUCCESS"
        for r in results
    )

    session_success = sum(
        r["session_status"] == "SUCCESS"
        for r in results
    )

    script_success = sum(
        r["script_status"] == "SUCCESS"
        for r in results
    )

    complete_success = sum(
        r["status"] == "SUCCESS"
        for r in results
    )

    logger.info("")
    logger.info("=" * 65)
    logger.info("SCRIPT 48 - FINAL SUMMARY")
    logger.info("=" * 65)

    logger.info("Proteins processed: %d", total)

    logger.info(
        "Cartoon PNG: %d/%d",
        cartoon_success,
        total,
    )

    logger.info(
        "Molecular Surface PNG: %d/%d",
        surface_success,
        total,
    )

    logger.info(
        "Surface + Cartoon PNG: %d/%d",
        combined_success,
        total,
    )

    logger.info(
        "PyMOL sessions: %d/%d",
        session_success,
        total,
    )

    logger.info(
        "PyMOL scripts: %d/%d",
        script_success,
        total,
    )

    logger.info(
        "Completely successful proteins: %d/%d",
        complete_success,
        total,
    )

    logger.info("")
    logger.info("AUDIT: %s", audit_path)
    logger.info("LOG: %s", LOG_FILE)

    logger.info("=" * 65)

    if complete_success != total:
        logger.warning(
            "Some proteins require attention. "
            "Check the audit CSV."
        )
    else:
        logger.info(
            "All proteins rendered successfully."
        )


if __name__ == "__main__":
    main()
