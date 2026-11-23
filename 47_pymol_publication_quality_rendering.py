
"""
TRANSCRIPTOME-X
SCRIPT 47 - PUBLICATION-QUALITY PROTEIN 3D RENDERING

Generates three PyMOL representations for each protein:
1. Cartoon / Ribbon
2. Molecular Surface
3. Transparent Surface + Cartoon

Also generates:
- High-resolution PNG images (300 DPI)
- Editable PyMOL sessions (.pse)
- PyMOL command scripts (.pml)
- Rendering audit CSV
- Execution log

Input:
Script 46 structural audit CSV.

Original structural files and Script 46 outputs are never modified.
"""

from __future__ import annotations

import csv
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path


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
    / "publication_3d_rendering"
)

CARTOON_DIR = OUTPUT_DIR / "cartoon"
SURFACE_DIR = OUTPUT_DIR / "molecular_surface"
COMBINED_DIR = OUTPUT_DIR / "surface_cartoon"

SESSION_DIR = OUTPUT_DIR / "sessions"
PML_DIR = OUTPUT_DIR / "pymol_scripts"
AUDIT_DIR = OUTPUT_DIR / "audit"

AUDIT_PATH = (
    AUDIT_DIR
    / "protein_publication_rendering_audit.csv"
)

LOG_PATH = AUDIT_DIR / "script47.log"

# Rendering dimensions.
IMAGE_WIDTH = 2400
IMAGE_HEIGHT = 2000
IMAGE_DPI = 300

# Molecular visualization.
SURFACE_TRANSPARENCY = 0.38

# PyMOL execution timeout.
PYMOL_TIMEOUT_SECONDS = 3600

# Optional explicit PyMOL executable.
# If empty, the script searches PATH and common locations.
PYMOL_EXECUTABLE = ""

# PyMOL executable candidates.
PYMOL_CANDIDATES = [
    r"C:\Program Files\PyMOL\PyMOL.exe",
    r"C:\Program Files\Schrodinger\PyMOL2\PyMOL.exe",
    r"C:\Program Files\Schrodinger\PyMOL2\bin\pymol.exe",
    r"C:\Program Files\Schrodinger\PyMOL\PyMOL.exe",
]

# Expected number of proteins.
EXPECTED_PROTEINS = 23


# ============================================================
# DIRECTORIES
# ============================================================

for directory in [
    CARTOON_DIR,
    SURFACE_DIR,
    COMBINED_DIR,
    SESSION_DIR,
    PML_DIR,
    AUDIT_DIR,
]:
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
# PYMOL DISCOVERY
# ============================================================

def find_pymol() -> str | None:

    if PYMOL_EXECUTABLE:

        executable = Path(PYMOL_EXECUTABLE)

        if executable.is_file():
            return str(executable)

        raise FileNotFoundError(
            f"Configured PyMOL executable not found: {executable}"
        )

    for name in [
        "pymol",
        "pymol.exe",
        "PyMOL.exe",
    ]:

        executable = shutil.which(name)

        if executable:
            return executable

    for candidate in PYMOL_CANDIDATES:

        path = Path(candidate)

        if path.is_file():
            return str(path)

    return None


def build_pymol_command(
    executable: str,
    script_path: Path,
) -> list[str]:

    suffix = Path(executable).suffix.lower()

    if suffix in {".bat", ".cmd"}:

        return [
            "cmd.exe",
            "/c",
            executable,
            "-cq",
            str(script_path),
        ]

    return [
        executable,
        "-cq",
        str(script_path),
    ]


# ============================================================
# INPUT AUDIT
# ============================================================

def load_structure_records() -> list[dict]:

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

        required_columns = {
            "gene",
            "uniprot_accession",
            "source",
            "structure_file",
            "status",
        }

        missing = required_columns - set(
            reader.fieldnames or []
        )

        if missing:

            raise ValueError(
                "Missing columns in Script 46 audit: "
                + ", ".join(sorted(missing))
            )

        for row in reader:

            gene = row["gene"].strip()
            structure_file = row["structure_file"].strip()

            if not gene or not structure_file:
                continue

            structure_path = Path(structure_file)

            if not structure_path.is_absolute():
                structure_path = PROJECT_DIR / structure_path

            if not structure_path.is_file():

                logger.warning(
                    "%s: structural file not found: %s",
                    gene,
                    structure_path,
                )

                continue

            if row["status"].strip() not in {
                "SUCCESS",
                "PARTIAL",
            }:

                logger.warning(
                    "%s: excluded because Script 46 status is %s",
                    gene,
                    row["status"],
                )

                continue

            records.append({
                "gene": gene,
                "accession": row["uniprot_accession"].strip(),
                "source": row["source"].strip(),
                "structure_file": structure_path,
            })

    if not records:

        raise RuntimeError(
            "No valid protein structures found in Script 46 audit."
        )

    genes = [
        item["gene"]
        for item in records
    ]

    if len(genes) != len(set(genes)):

        raise ValueError(
            "Duplicate gene entries found in Script 46 audit."
        )

    logger.info(
        "Valid structures found: %s",
        len(records),
    )

    return records


# ============================================================
# PYMOL STRING UTILITIES
# ============================================================

def pymol_quote(value: str | Path) -> str:

    normalized = str(value).replace("\\", "/")
    normalized = normalized.replace('"', '\\"')

    return f'"{normalized}"'


# ============================================================
# PYMOL GLOBAL SETTINGS
# ============================================================

def build_global_settings() -> list[str]:

    return [
        "reinitialize",
        "set retain_order, 1",
        "set auto_zoom, 0",
        "set depth_cue, 0",
        "set ray_opaque_background, 1",
        "set ray_shadows, 1",
        "set ray_trace_fog, 0",
        "set antialias, 2",
        "set ambient, 0.35",
        "set direct, 0.75",
        "set specular, 0.25",
        "set shininess, 35",
        "set reflect, 0.08",
        "set cartoon_sampling, 14",
        "set cartoon_smooth_loops, 1",
        "set cartoon_fancy_helices, 1",
        "set cartoon_flat_sheets, 1",
        "set cartoon_transparency, 0",
        "set surface_quality, 2",
        "set transparency_mode, 2",
        "bg_color white",
        "set ray_trace_mode, 0",
    ]


# ============================================================
# STRUCTURE LOADING
# ============================================================

def build_load_commands(
    record: dict,
    object_name: str,
) -> list[str]:

    path = pymol_quote(
        record["structure_file"]
    )

    return [
        f"load {path}, {object_name}",
        f"remove solvent, {object_name}",
        f"dss {object_name}",
        f"hide everything, {object_name}",
    ]


# ============================================================
# CARTOON COLORING
# ============================================================

def build_cartoon_commands(
    object_name: str,
) -> list[str]:

    selection = object_name

    return [
        f"show cartoon, {selection}",
        f"color gray70, {selection}",
        f"color salmon, {selection} and ss H",
        f"color skyblue, {selection} and ss S",
        f"color gray70, {selection} and ss L",
        f"set cartoon_transparency, 0, {selection}",
        f"set cartoon_flat_sheets, 1",
    ]


# ============================================================
# VIEW AND RENDER
# ============================================================

def build_render_commands(
    object_name: str,
    output_path: Path,
) -> list[str]:

    return [
        f"orient {object_name}",
        f"zoom {object_name}, 4",
        f"ray {IMAGE_WIDTH}, {IMAGE_HEIGHT}",
        (
            f"png {pymol_quote(output_path)}, "
            f"dpi={IMAGE_DPI}"
        ),
    ]


# ============================================================
# CARTOON RENDERING
# ============================================================

def build_cartoon_render(
    object_name: str,
    output_path: Path,
) -> list[str]:

    commands = [
        f"hide everything, {object_name}",
        f"show cartoon, {object_name}",
    ]

    commands.extend(
        build_cartoon_commands(object_name)
    )

    commands.extend(
        build_render_commands(
            object_name,
            output_path,
        )
    )

    return commands


# ============================================================
# MOLECULAR SURFACE RENDERING
# ============================================================

def build_surface_render(
    object_name: str,
    output_path: Path,
) -> list[str]:

    return [
        f"hide everything, {object_name}",
        f"show surface, {object_name}",
        f"color gray80, {object_name}",
        f"set surface_quality, 2",
        f"set transparency, 0, {object_name}",
        f"set specular, 0.35",
        f"set shininess, 40",
        *build_render_commands(
            object_name,
            output_path,
        ),
    ]


# ============================================================
# SURFACE + CARTOON RENDERING
# ============================================================

def build_combined_render(
    object_name: str,
    output_path: Path,
) -> list[str]:

    commands = [
        f"hide everything, {object_name}",

        # Molecular surface.
        f"show surface, {object_name}",
        f"color lightblue, {object_name}",
        (
            f"set transparency, "
            f"{SURFACE_TRANSPARENCY}, {object_name}"
        ),

        # Internal protein structure.
        f"show cartoon, {object_name}",
        f"set cartoon_transparency, 0, {object_name}",
    ]

    commands.extend(
        build_cartoon_commands(object_name)
    )

    commands.extend(
        build_render_commands(
            object_name,
            output_path,
        )
    )

    return commands


# ============================================================
# PYMOL SCRIPT GENERATION
# ============================================================

def generate_pymol_script(
    records: list[dict],
) -> Path:

    script_path = (
        PML_DIR
        / "script47_render_all_proteins.pml"
    )

    commands = []

    commands.extend(
        build_global_settings()
    )

    for index, record in enumerate(
        records,
        start=1,
    ):

        gene = record["gene"]

        object_name = f"protein_{index:02d}"

        cartoon_path = (
            CARTOON_DIR
            / f"{gene}_cartoon.png"
        )

        surface_path = (
            SURFACE_DIR
            / f"{gene}_surface.png"
        )

        combined_path = (
            COMBINED_DIR
            / f"{gene}_surface_cartoon.png"
        )

        session_path = (
            SESSION_DIR
            / f"{gene}.pse"
        )

        commands.append(
            f"\n# ========================================"
        )

        commands.append(
            f"# {gene} | {record['accession']}"
        )

        commands.append(
            f"# Source: {record['source']}"
        )

        commands.append(
            f"# ========================================"
        )

        commands.extend(
            build_load_commands(
                record,
                object_name,
            )
        )

        # 1. Cartoon / Ribbon.
        commands.extend(
            build_cartoon_render(
                object_name,
                cartoon_path,
            )
        )

        # 2. Molecular surface.
        commands.extend(
            build_surface_render(
                object_name,
                surface_path,
            )
        )

        # 3. Surface + Cartoon.
        commands.extend(
            build_combined_render(
                object_name,
                combined_path,
            )
        )

        # Save editable PyMOL session in combined view.
        commands.append(
            f"save {pymol_quote(session_path)}"
        )

        commands.append(
            f"delete {object_name}"
        )

    commands.extend([
        "bg_color white",
        "quit",
    ])

    script_path.write_text(
        "\n".join(commands),
        encoding="utf-8",
    )

    return script_path


# ============================================================
# AUDIT
# ============================================================

def build_audit_records(
    records: list[dict],
) -> list[dict]:

    audit_records = []

    for record in records:

        gene = record["gene"]

        paths = {
            "cartoon": (
                CARTOON_DIR
                / f"{gene}_cartoon.png"
            ),
            "surface": (
                SURFACE_DIR
                / f"{gene}_surface.png"
            ),
            "combined": (
                COMBINED_DIR
                / f"{gene}_surface_cartoon.png"
            ),
            "session": (
                SESSION_DIR
                / f"{gene}.pse"
            ),
        }

        audit_records.append({
            "gene": gene,
            "uniprot_accession": record["accession"],
            "source": record["source"],
            "structure_file": str(record["structure_file"]),
            "cartoon_png": str(paths["cartoon"]),
            "surface_png": str(paths["surface"]),
            "combined_png": str(paths["combined"]),
            "pymol_session": str(paths["session"]),
            "cartoon_status": "PENDING",
            "surface_status": "PENDING",
            "combined_status": "PENDING",
            "session_status": "PENDING",
            "status": "PENDING",
            "error": "",
        })

    return audit_records


def validate_output(path: Path) -> bool:

    return (
        path.is_file()
        and path.stat().st_size > 1000
    )


def update_audit_status(
    audit_records: list[dict],
) -> None:

    for record in audit_records:

        gene = record["gene"]

        outputs = {
            "cartoon_status": (
                CARTOON_DIR
                / f"{gene}_cartoon.png"
            ),
            "surface_status": (
                SURFACE_DIR
                / f"{gene}_surface.png"
            ),
            "combined_status": (
                COMBINED_DIR
                / f"{gene}_surface_cartoon.png"
            ),
            "session_status": (
                SESSION_DIR
                / f"{gene}.pse"
            ),
        }

        for field, path in outputs.items():

            record[field] = (
                "SUCCESS"
                if validate_output(path)
                else "FAILED"
            )

        statuses = [
            record["cartoon_status"],
            record["surface_status"],
            record["combined_status"],
            record["session_status"],
        ]

        if all(status == "SUCCESS" for status in statuses):
            record["status"] = "SUCCESS"

        elif any(status == "SUCCESS" for status in statuses):
            record["status"] = "PARTIAL"

        else:
            record["status"] = "FAILED"


def save_audit(
    audit_records: list[dict],
) -> None:

    if not audit_records:
        return

    with AUDIT_PATH.open(
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


# ============================================================
# MAIN
# ============================================================

def main():

    logger.info("=" * 75)
    logger.info("TRANSCRIPTOME-X - SCRIPT 47")
    logger.info("PUBLICATION-QUALITY PROTEIN 3D RENDERING")
    logger.info("=" * 75)

    records = load_structure_records()

    if len(records) != EXPECTED_PROTEINS:

        logger.warning(
            "Expected %s proteins, found %s.",
            EXPECTED_PROTEINS,
            len(records),
        )

    logger.info(
        "Structures to render: %s",
        len(records),
    )

    logger.info(
        "Output directory: %s",
        OUTPUT_DIR,
    )

    # Prepare audit.
    audit_records = build_audit_records(records)

    # Generate the PyMOL command script.
    pml_path = generate_pymol_script(records)

    logger.info(
        "PyMOL script created: %s",
        pml_path,
    )

    # Find PyMOL.
    executable = find_pymol()

    if executable is None:

        for record in audit_records:

            record["status"] = "PYMOL_NOT_FOUND"
            record["error"] = (
                "PyMOL executable not detected. "
                "PML script generated but rendering not executed."
            )

        save_audit(audit_records)

        print("\n" + "=" * 75)
        print("PYMOL NOT FOUND")
        print("=" * 75)

        print(
            "\nThe rendering script has been generated, "
            "but PyMOL is not available."
        )

        print("\nGenerated PML:")
        print(pml_path)

        print("\nInstall PyMOL or configure PYMOL_EXECUTABLE.")
        print("Then execute this script again.")

        print("\nAudit:")
        print(AUDIT_PATH)

        print("=" * 75)

        return

    logger.info(
        "PyMOL executable: %s",
        executable,
    )

    command = build_pymol_command(
        executable,
        pml_path,
    )

    logger.info(
        "Starting PyMOL rendering process..."
    )

    try:

        result = subprocess.run(
            command,
            cwd=str(PROJECT_DIR),
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=PYMOL_TIMEOUT_SECONDS,
            check=False,
        )

        pymol_log_path = (
            AUDIT_DIR
            / "pymol_execution.log"
        )

        pymol_log_path.write_text(
            "STDOUT\n"
            + (result.stdout or "")
            + "\n\nSTDERR\n"
            + (result.stderr or ""),
            encoding="utf-8",
        )

        if result.returncode != 0:

            logger.error(
                "PyMOL returned exit code %s.",
                result.returncode,
            )

            logger.error(
                "Check: %s",
                pymol_log_path,
            )

            for record in audit_records:
                record["error"] = (
                    f"PyMOL exit code: {result.returncode}"
                )

        else:

            logger.info(
                "PyMOL process completed."
            )

    except subprocess.TimeoutExpired:

        logger.exception(
            "PyMOL rendering exceeded the timeout."
        )

        for record in audit_records:
            record["error"] = "PYMOL_TIMEOUT"

    except Exception as exc:

        logger.exception(
            "PyMOL execution failed."
        )

        for record in audit_records:
            record["error"] = str(exc)

    # Validate every output.
    update_audit_status(audit_records)

    save_audit(audit_records)

    # ========================================================
    # FINAL SUMMARY
    # ========================================================

    cartoon_success = sum(
        item["cartoon_status"] == "SUCCESS"
        for item in audit_records
    )

    surface_success = sum(
        item["surface_status"] == "SUCCESS"
        for item in audit_records
    )

    combined_success = sum(
        item["combined_status"] == "SUCCESS"
        for item in audit_records
    )

    session_success = sum(
        item["session_status"] == "SUCCESS"
        for item in audit_records
    )

    complete = sum(
        item["status"] == "SUCCESS"
        for item in audit_records
    )

    print("\n" + "=" * 75)
    print("SCRIPT 47 - FINAL SUMMARY")
    print("=" * 75)

    print(f"Proteins processed: {len(records)}")
    print(f"Cartoon PNG: {cartoon_success}/{len(records)}")
    print(f"Molecular Surface PNG: {surface_success}/{len(records)}")
    print(f"Surface + Cartoon PNG: {combined_success}/{len(records)}")
    print(f"PyMOL sessions: {session_success}/{len(records)}")

    print(
        f"\nCompletely successful proteins: "
        f"{complete}/{len(records)}"
    )

    print("\nRESULTS BY PROTEIN")

    for item in audit_records:

        print(
            f"{item['gene']:10s} | "
            f"{item['status']:12s} | "
            f"Cartoon: {item['cartoon_status']:7s} | "
            f"Surface: {item['surface_status']:7s} | "
            f"Combined: {item['combined_status']:7s}"
        )

    print("\nOUTPUT DIRECTORIES")
    print(f"Cartoon:  {CARTOON_DIR}")
    print(f"Surface:  {SURFACE_DIR}")
    print(f"Combined: {COMBINED_DIR}")
    print(f"Sessions: {SESSION_DIR}")
    print(f"Audit:    {AUDIT_PATH}")

    print("=" * 75)


if __name__ == "__main__":
    main()
