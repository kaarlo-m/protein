from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List

import MDAnalysis as mda

BASE_DATA_DIR = Path("data")
DEFAULT_FRAME_STRIDE = 10


def iter_analysis_dirs(base_dir: Path) -> Iterable[Path]:
    """Yield analysis directories located directly under the base data directory."""
    for child in sorted(base_dir.iterdir()):
        if child.is_dir() and child.name.endswith("_analysis"):
            yield child


def find_structure_file(analysis_dir: Path) -> Path:
    """Return the first PDB file inside an analysis directory."""
    pdb_files = sorted(analysis_dir.glob("*.pdb"))
    if not pdb_files:
        raise FileNotFoundError(f"No structure (.pdb) file found in {analysis_dir}")
    return pdb_files[0]


def write_multimodel(structure: Path, trajectory: Path, stride: int) -> Path:
    """Generate a multi-model PDB trajectory sampled with the provided stride."""
    universe = mda.Universe(str(structure), str(trajectory))
    atoms = universe.select_atoms("protein")
    output_path = trajectory.with_name(f"{trajectory.stem}_multimodel.pdb")

    with mda.Writer(str(output_path), atoms.n_atoms) as writer:
        for ts in universe.trajectory[::stride]:
            writer.write(atoms)

    return output_path


def process_analysis_dir(analysis_dir: Path, stride: int) -> List[Path]:
    structure = find_structure_file(analysis_dir)
    trajectories = sorted(analysis_dir.glob("*.xtc"))
    created_files = []

    if not trajectories:
        print(f"[WARN] No trajectories found in {analysis_dir}")
        return created_files

    # Prefer the R1 trajectory if present; otherwise fall back to the first available.
    first_traj = next((t for t in trajectories if "_R1" in t.stem), trajectories[0])
    output_path = first_traj.with_name(f"{first_traj.stem}_multimodel.pdb")

    if output_path.exists():
        print(f"[SKIP] Multi-model already exists: {output_path}")
        return created_files

    print(f"[INFO] Writing multi-model PDB for {first_traj}")
    output = write_multimodel(structure, first_traj, stride)
    created_files.append(output)
    print(f"[OK]   Saved to {output}")

    return created_files


def main(base_dir: Path, stride: int) -> None:
    if not base_dir.exists():
        raise FileNotFoundError(f"Base data directory not found: {base_dir}")

    for analysis_dir in iter_analysis_dirs(base_dir):
        print(f"[INFO] Processing {analysis_dir}")
        process_analysis_dir(analysis_dir, stride)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate multi-model PDB files for all analysis trajectories."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=BASE_DATA_DIR,
        help="Path to the data directory containing *_analysis folders.",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=DEFAULT_FRAME_STRIDE,
        help="Frame stride when sampling trajectories.",
    )

    args = parser.parse_args()
    main(args.data_dir, args.stride)
