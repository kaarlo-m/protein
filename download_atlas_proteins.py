from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path
from zipfile import ZipFile, BadZipFile

DEFAULT_PROTEINS_FILE = Path("proteins_list.csv")
DEFAULT_DATA_DIR = Path("data")
API_TEMPLATE = "https://www.dsimb.inserm.fr/ATLAS/api/ATLAS/analysis/{pdb_id}"


def read_pdb_ids(csv_path: Path) -> list[str]:
    """Read PDB identifiers from proteins_list.csv."""
    if not csv_path.exists():
        raise FileNotFoundError(f"Proteins file not found: {csv_path}")

    pdb_ids: list[str] = []
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            pdb = row.get("PDB", "").strip()
            if pdb:
                pdb_ids.append(pdb)
    return pdb_ids


def download_archive(pdb_id: str, dest_zip: Path) -> None:
    """Download a single PDB archive via curl."""
    url = API_TEMPLATE.format(pdb_id=pdb_id)
    cmd = [
        "curl",
        "-X",
        "GET",
        url,
        "-H",
        "accept: application/octet-stream",
        "--output",
        str(dest_zip),
        "--fail",
        "--location",
    ]
    print(f"[INFO] Downloading {pdb_id} -> {dest_zip}")
    subprocess.run(cmd, check=True)


def extract_archive(zip_path: Path, target_dir: Path) -> None:
    """Extract the archive into target_dir."""
    target_dir.mkdir(parents=True, exist_ok=True)
    try:
        with ZipFile(zip_path, "r") as zf:
            zf.extractall(target_dir)
    except BadZipFile as exc:
        raise RuntimeError(f"Corrupt zip file {zip_path}") from exc


def process_protein(pdb_id: str, data_dir: Path, force: bool) -> None:
    """Download and extract a protein archive."""
    analysis_dir = data_dir / f"{pdb_id}_analysis"
    zip_path = analysis_dir.with_suffix(".zip")

    if analysis_dir.exists() and not force:
        print(f"[SKIP] Analysis folder already present for {pdb_id}: {analysis_dir}")
        return

    analysis_dir.mkdir(parents=True, exist_ok=True)

    try:
        download_archive(pdb_id, zip_path)
    except subprocess.CalledProcessError as exc:
        print(f"[ERROR] Failed to download {pdb_id}: {exc}")
        return

    try:
        extract_archive(zip_path, analysis_dir)
    except RuntimeError as exc:
        print(f"[ERROR] Extract failed for {pdb_id}: {exc}")
        return
    finally:
        if zip_path.exists():
            zip_path.unlink()

    print(f"[OK]   Extracted {pdb_id} to {analysis_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Download protein analyses from ATLAS and extract to data/[PDB]_analysis."
    )
    parser.add_argument(
        "--proteins-file",
        type=Path,
        default=DEFAULT_PROTEINS_FILE,
        help="CSV file containing a PDB column (default: proteins_list.csv).",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help="Destination data directory (default: data).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-download and overwrite even if the analysis directory exists.",
    )
    args = parser.parse_args()

    try:
        pdb_ids = read_pdb_ids(args.proteins_file)
    except FileNotFoundError as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)

    if not pdb_ids:
        print("[WARN] No PDB IDs found to process.")
        return

    args.data_dir.mkdir(parents=True, exist_ok=True)

    for pdb_id in pdb_ids:
        process_protein(pdb_id, args.data_dir, args.force)


if __name__ == "__main__":
    main()
