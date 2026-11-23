
"""
TRANSCRIPTOME-X
SCRIPT 45 - RECOVERY OF MISSING PROTEIN STRUCTURES

Obiettivo:
- Analizzare le 23 proteine del progetto.
- Preservare tutte le strutture AlphaFold già scaricate.
- Recuperare le strutture mancanti da AlphaFold DB.
- Cercare strutture sperimentali nel PDB.
- Scaricare la sequenza FASTA da UniProt.
- Preparare gli input FASTA per ColabFold quando necessario.
- Produrre un audit CSV e un report finale.

Le strutture sperimentali PDB e quelle predette non vengono mai
considerate equivalenti e sono salvate in directory separate.
"""

from __future__ import annotations

import csv
import json
import logging
import re
import sys
import time
from pathlib import Path
from typing import Any

import requests


# ============================================================
# CONFIGURATION
# ============================================================

PROJECT_DIR = Path(
    r"C:\Users\luank\Dropbox\PC\Desktop\PYTHON\TRANSCRIPTOME_X"
)

BASE_DIR = PROJECT_DIR / "results" / "gene_protein_illustrations"

AUDIT_FILE = BASE_DIR / "audit" / "gene_protein_illustration_audit.csv"

EXISTING_AF_DIR = BASE_DIR / "alphafold_3d"

RECOVERY_DIR = BASE_DIR / "structure_recovery"

ALPHAFOLD_DIR = RECOVERY_DIR / "alphafold"
PDB_DIR = RECOVERY_DIR / "experimental_pdb"
FASTA_DIR = RECOVERY_DIR / "fasta"
COLABFOLD_DIR = RECOVERY_DIR / "colabfold_input"
REPORT_DIR = RECOVERY_DIR / "audit"

UNIPROT_BASE = "https://rest.uniprot.org"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"
RCSB_SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2/query"

REQUEST_TIMEOUT = 60
RETRY_COUNT = 3
REQUEST_DELAY = 1.0

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


# ============================================================
# LOGGING
# ============================================================

RECOVERY_DIR.mkdir(parents=True, exist_ok=True)

for directory in [
    ALPHAFOLD_DIR,
    PDB_DIR,
    FASTA_DIR,
    COLABFOLD_DIR,
    REPORT_DIR,
]:
    directory.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
    handlers=[
        logging.FileHandler(
            REPORT_DIR / "script45_recovery.log",
            encoding="utf-8",
        ),
        logging.StreamHandler(sys.stdout),
    ],
)

logger = logging.getLogger("TRANSCRIPTOME-X")


# ============================================================
# HTTP SESSION
# ============================================================

SESSION = requests.Session()

SESSION.headers.update({
    "User-Agent": "TRANSCRIPTOME-X/1.0 (protein structure recovery)",
    "Accept": "application/json, text/plain, */*",
})


def request_with_retry(
    method: str,
    url: str,
    **kwargs,
) -> requests.Response | None:

    for attempt in range(1, RETRY_COUNT + 1):

        try:
            response = SESSION.request(
                method,
                url,
                timeout=REQUEST_TIMEOUT,
                **kwargs,
            )

            if response.status_code == 404:
                return None

            response.raise_for_status()
            return response

        except requests.RequestException as exc:

            logger.warning(
                "Request failed (%s/%s): %s | %s",
                attempt,
                RETRY_COUNT,
                url,
                exc,
            )

            if attempt < RETRY_COUNT:
                time.sleep(2 * attempt)

    return None


# ============================================================
# FILE UTILITIES
# ============================================================

def safe_write_text(path: Path, content: str) -> None:

    path.parent.mkdir(parents=True, exist_ok=True)

    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(content, encoding="utf-8")
    temporary.replace(path)


def download_file(url: str, destination: Path) -> bool:

    if destination.exists() and destination.stat().st_size > 100:
        logger.info("Already present: %s", destination.name)
        return True

    response = request_with_retry("GET", url)

    if response is None:
        return False

    if len(response.content) < 100:
        return False

    destination.parent.mkdir(parents=True, exist_ok=True)

    temporary = destination.with_suffix(
        destination.suffix + ".tmp"
    )

    temporary.write_bytes(response.content)
    temporary.replace(destination)

    return True


def load_original_audit() -> dict[str, dict[str, str]]:

    if not AUDIT_FILE.exists():
        logger.warning("Original audit not found: %s", AUDIT_FILE)
        return {}

    records = {}

    with AUDIT_FILE.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        for row in csv.DictReader(handle):
            records[row["gene"]] = row

    return records


def existing_structure_available(gene: str) -> bool:

    if not EXISTING_AF_DIR.exists():
        return False

    patterns = [
        f"*{gene}*.cif",
        f"*{gene}*.pdb",
        f"*{gene}*.html",
    ]

    for pattern in patterns:
        if list(EXISTING_AF_DIR.rglob(pattern)):
            return True

    return False


# ============================================================
# UNIPROT FASTA
# ============================================================

def retrieve_uniprot_fasta(
    gene: str,
    accession: str,
) -> tuple[str | None, str]:

    url = f"{UNIPROT_BASE}/uniprotkb/{accession}.fasta"

    response = request_with_retry("GET", url)

    if response is None:
        return None, "UNAVAILABLE"

    fasta = response.text.strip()

    if not fasta.startswith(">"):
        return None, "INVALID_FASTA"

    destination = FASTA_DIR / f"{gene}_{accession}.fasta"

    safe_write_text(destination, fasta + "\n")

    sequence = "".join(
        line.strip()
        for line in fasta.splitlines()
        if not line.startswith(">")
    )

    if not sequence:
        return None, "EMPTY_SEQUENCE"

    logger.info(
        "%s: UniProt FASTA retrieved (%s aa)",
        gene,
        len(sequence),
    )

    return fasta, "SUCCESS"


# ============================================================
# ALPHAFOLD API
# ============================================================

def extract_structure_urls(data: Any) -> list[dict[str, str]]:
    """
    Estrae URL di strutture dai record API senza dipendere
    da una singola versione dello schema JSON.
    """

    results = []

    def walk(obj: Any):

        if isinstance(obj, dict):

            normalized = {
                str(key).lower(): value
                for key, value in obj.items()
            }

            for key, value in normalized.items():

                if (
                    isinstance(value, str)
                    and value.startswith("https://")
                    and any(
                        token in key
                        for token in [
                            "cifurl",
                            "pdburl",
                            "bcifurl",
                            "modelurl",
                        ]
                    )
                ):
                    results.append({
                        "format": (
                            "cif" if "cif" in key
                            else "pdb" if "pdb" in key
                            else "bcif"
                        ),
                        "url": value,
                    })

            for value in obj.values():
                walk(value)

        elif isinstance(obj, list):

            for item in obj:
                walk(item)

    walk(data)

    unique = {}

    for item in results:
        unique[item["url"]] = item

    return list(unique.values())


def retrieve_alphafold(
    gene: str,
    accession: str,
) -> tuple[int, str]:

    url = f"{ALPHAFOLD_API}/{accession}"

    response = request_with_retry("GET", url)

    if response is None:
        return 0, "UNAVAILABLE"

    try:
        data = response.json()
    except ValueError:
        return 0, "INVALID_API_RESPONSE"

    structure_urls = extract_structure_urls(data)

    # Prefer mmCIF, then PDB. Binary CIF is not downloaded
    # because the current workflow uses ordinary CIF/PDB files.
    structure_urls = [
        item for item in structure_urls
        if item["format"] in ("cif", "pdb")
    ]

    structure_urls.sort(
        key=lambda item: 0 if item["format"] == "cif" else 1
    )

    downloaded = 0
    attempted_urls = set()

    for index, item in enumerate(structure_urls, start=1):

        structure_url = item["url"]

        if structure_url in attempted_urls:
            continue

        attempted_urls.add(structure_url)

        extension = item["format"]

        filename = (
            f"{gene}_{accession}_AlphaFold_{index}.{extension}"
        )

        destination = ALPHAFOLD_DIR / gene / filename

        if download_file(structure_url, destination):
            downloaded += 1
            logger.info(
                "%s: AlphaFold structure downloaded: %s",
                gene,
                filename,
            )

        time.sleep(REQUEST_DELAY)

    if downloaded:
        return downloaded, "DOWNLOADED"

    if not structure_urls:
        return 0, "NO_STRUCTURE_URL_RETURNED"

    return 0, "DOWNLOAD_FAILED"


# ============================================================
# RCSB PDB SEARCH
# ============================================================

def search_rcsb(accession: str) -> list[str]:

    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": (
                    "rcsb_polymer_entity_container_identifiers."
                    "reference_sequence_identifiers.database_accession"
                ),
                "operator": "exact_match",
                "value": accession,
            },
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": 100,
            },
            "results_content_type": ["experimental"],
        },
    }

    response = request_with_retry(
        "POST",
        RCSB_SEARCH_API,
        json=query,
        headers={"Content-Type": "application/json"},
    )

    if response is None:
        return []

    try:
        data = response.json()
    except ValueError:
        return []

    identifiers = []

    for item in data.get("result_set", []):
        identifier = item.get("identifier")

        if identifier:
            identifiers.append(identifier)

    return sorted(set(identifiers))


def download_rcsb_structures(
    gene: str,
    accession: str,
) -> tuple[int, list[str]]:

    identifiers = search_rcsb(accession)

    if not identifiers:
        return 0, []

    downloaded = 0
    saved_ids = []

    for pdb_id in identifiers:

        pdb_id = pdb_id.upper()

        url = (
            f"https://files.rcsb.org/download/{pdb_id}.cif"
        )

        destination = (
            PDB_DIR / gene / f"{gene}_{pdb_id}_experimental.cif"
        )

        if download_file(url, destination):

            downloaded += 1
            saved_ids.append(pdb_id)

            logger.info(
                "%s: experimental PDB structure retrieved: %s",
                gene,
                pdb_id,
            )

        time.sleep(REQUEST_DELAY)

    return downloaded, saved_ids


# ============================================================
# COLABFOLD INPUT
# ============================================================

def prepare_colabfold_input(
    gene: str,
    fasta: str | None,
) -> str:

    if not fasta:
        return "NOT_PREPARED_NO_FASTA"

    # ColabFold accepts FASTA input directly.
    # This file is prepared for remote or local execution.
    destination = COLABFOLD_DIR / f"{gene}_colabfold.fasta"

    safe_write_text(destination, fasta + "\n")

    return "READY"


# ============================================================
# MAIN RECOVERY WORKFLOW
# ============================================================

def main():

    original_audit = load_original_audit()

    results = []

    logger.info("=" * 65)
    logger.info("TRANSCRIPTOME-X - SCRIPT 45")
    logger.info("PROTEIN STRUCTURE RECOVERY")
    logger.info("Total genes: %s", len(GENES))
    logger.info("=" * 65)

    for position, (gene, accession) in enumerate(
        GENES,
        start=1,
    ):

        logger.info("")
        logger.info(
            "[%s/%s] Processing %s (%s)",
            position,
            len(GENES),
            gene,
            accession,
        )

        original = original_audit.get(gene, {})

        original_status = original.get(
            "structure_status",
            "UNKNOWN",
        )

        original_af_status = original.get(
            "alphafold_status",
            "UNKNOWN",
        )

        # Preserve all existing structures.
        if (
            original_status == "SUCCESS"
            and original_af_status == "DOWNLOADED"
        ):

            logger.info(
                "%s: original AlphaFold structure already available. "
                "Skipping recovery.",
                gene,
            )

            results.append({
                "gene": gene,
                "uniprot_accession": accession,
                "original_structure_status": original_status,
                "original_alphafold_status": original_af_status,
                "existing_structure_detected": (
                    existing_structure_available(gene)
                ),
                "alphafold_recovery": "NOT_REQUIRED",
                "alphafold_models_downloaded": 0,
                "pdb_structures_downloaded": 0,
                "pdb_ids": "",
                "fasta_status": "NOT_REQUIRED",
                "fasta_length": "",
                "colabfold_status": "NOT_REQUIRED",
                "final_status": "ORIGINAL_MODEL_AVAILABLE",
            })

            continue

        # ----------------------------------------------------
        # 1. AlphaFold
        # ----------------------------------------------------

        af_count, af_status = retrieve_alphafold(
            gene,
            accession,
        )

        # ----------------------------------------------------
        # 2. UniProt FASTA
        # ----------------------------------------------------

        fasta, fasta_status = retrieve_uniprot_fasta(
            gene,
            accession,
        )

        sequence_length = ""

        if fasta:
            sequence_length = len(
                "".join(
                    line.strip()
                    for line in fasta.splitlines()
                    if not line.startswith(">")
                )
            )

        # ----------------------------------------------------
        # 3. Experimental PDB
        # ----------------------------------------------------

        pdb_count, pdb_ids = download_rcsb_structures(
            gene,
            accession,
        )

        # ----------------------------------------------------
        # 4. ColabFold FASTA
        # ----------------------------------------------------

        colabfold_status = "NOT_REQUIRED"

        if af_count == 0:

            colabfold_status = prepare_colabfold_input(
                gene,
                fasta,
            )

        # ----------------------------------------------------
        # 5. Final classification
        # ----------------------------------------------------

        if af_count > 0:
            final_status = "ALPHAFOLD_MODEL_RECOVERED"

        elif pdb_count > 0:
            final_status = "EXPERIMENTAL_PDB_AVAILABLE"

        elif fasta:
            final_status = "FASTA_READY_FOR_COLABFOLD"

        else:
            final_status = "NO_STRUCTURE_OR_SEQUENCE_AVAILABLE"

        results.append({
            "gene": gene,
            "uniprot_accession": accession,
            "original_structure_status": original_status,
            "original_alphafold_status": original_af_status,
            "existing_structure_detected": (
                existing_structure_available(gene)
            ),
            "alphafold_recovery": af_status,
            "alphafold_models_downloaded": af_count,
            "pdb_structures_downloaded": pdb_count,
            "pdb_ids": ";".join(pdb_ids),
            "fasta_status": fasta_status,
            "fasta_length": sequence_length,
            "colabfold_status": colabfold_status,
            "final_status": final_status,
        })

        logger.info(
            "%s FINAL STATUS: %s",
            gene,
            final_status,
        )

    # ========================================================
    # SAVE REPORTS
    # ========================================================

    report_csv = REPORT_DIR / "structure_recovery_report.csv"

    if results:

        fieldnames = list(results[0].keys())

        with report_csv.open(
            "w",
            encoding="utf-8-sig",
            newline="",
        ) as handle:

            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
            )

            writer.writeheader()
            writer.writerows(results)

    summary = {
        "project": "TRANSCRIPTOME-X",
        "script": "45_recover_missing_protein_structures.py",
        "total_genes": len(GENES),
        "original_models_preserved": sum(
            1 for item in results
            if item["final_status"] == "ORIGINAL_MODEL_AVAILABLE"
        ),
        "new_alphafold_models": sum(
            item["alphafold_models_downloaded"]
            for item in results
        ),
        "experimental_pdb_structures": sum(
            item["pdb_structures_downloaded"]
            for item in results
        ),
        "fasta_sequences_available": sum(
            1 for item in results
            if item["fasta_status"] == "SUCCESS"
        ),
        "colabfold_inputs_ready": sum(
            1 for item in results
            if item["colabfold_status"] == "READY"
        ),
        "results": results,
    }

    summary_json = REPORT_DIR / "structure_recovery_summary.json"

    summary_json.write_text(
        json.dumps(
            summary,
            indent=4,
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    # ========================================================
    # FINAL CONSOLE REPORT
    # ========================================================

    print("\n" + "=" * 65)
    print("TRANSCRIPTOME-X - SCRIPT 45 COMPLETED")
    print("=" * 65)

    print(f"Total genes: {len(GENES)}")

    print(
        "Original structures preserved:",
        summary["original_models_preserved"],
    )

    print(
        "New AlphaFold models:",
        summary["new_alphafold_models"],
    )

    print(
        "Experimental PDB structures:",
        summary["experimental_pdb_structures"],
    )

    print(
        "FASTA sequences available:",
        summary["fasta_sequences_available"],
    )

    print(
        "ColabFold inputs ready:",
        summary["colabfold_inputs_ready"],
    )

    print("\nFINAL STATUS BY GENE")

    for item in results:

        print(
            f"{item['gene']:10s} | "
            f"{item['final_status']}"
        )

    print("\nOUTPUT DIRECTORIES")
    print(f"AlphaFold: {ALPHAFOLD_DIR}")
    print(f"Experimental PDB: {PDB_DIR}")
    print(f"FASTA: {FASTA_DIR}")
    print(f"ColabFold: {COLABFOLD_DIR}")
    print(f"Audit: {REPORT_DIR}")

    print("\nReports:")
    print(report_csv)
    print(summary_json)

    print("=" * 65)


if __name__ == "__main__":
    main()
