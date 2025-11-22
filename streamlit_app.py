from pathlib import Path

import csv
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components

DATA_DIR = Path("data")
ATLAS_INFO_PATH = Path("ATLAS_info.tsv")
VIEW_WIDTH = 1000
VIEW_HEIGHT = 1000

st.set_page_config(
    layout="wide"
)


def discover_multimodel_files(base_dir: Path = DATA_DIR) -> list[Path]:
    """Return all multi-model PDB files located under the data directory."""
    if not base_dir.exists():
        return []
    return sorted(base_dir.rglob("*_multimodel.pdb"))


@st.cache_data(show_spinner=False)
def load_multimodel_pdb(pdb_path: str) -> str:
    """Read and cache the contents of a multi-model PDB file."""
    return Path(pdb_path).read_text()


@st.cache_data(show_spinner=False)
def load_atlas_metadata(tsv_path: Path = ATLAS_INFO_PATH) -> dict[str, dict[str, str]]:
    """Parse ATLAS_info.tsv for quick protein lookups by PDB id."""
    if not tsv_path.exists():
        return {}

    records: dict[str, dict[str, str]] = {}
    with tsv_path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            pdb_id = row.get("PDB", "").strip()
            if not pdb_id:
                continue
            records[pdb_id] = {
                "name": pdb_id,
                "title": row.get("protein_name", "").strip(),
                "organism": row.get("organism", "").strip(),
                "sequence": row.get("sequence", "").strip(),
                "uniprot" : row.get("UniProt_entry", "").strip(),
            }
    return records


def render_mol(pdb_data: str) -> None:
    view = py3Dmol.view(width=VIEW_WIDTH, height=VIEW_HEIGHT)
    view.addModelsAsFrames(pdb_data, "pdb")
    view.setStyle({"cartoon": {"color": "spectrum"}})
    view.animate({"interval": 75})
    view.setBackgroundColor("black")
    view.zoomTo()
    view.spin(True)

    viewer_markup = view._make_html()
    framed_html = f"""
    <style>
        body {{
            margin: 0;
            background-color: transparent;
        }}
        .viewer-frame {{
            border: 2px solid #fbff00ff;
            border-radius: 12px;
            padding: 12px;
            background-color: #050505;
            box-shadow: 0 0 20px rgba(76, 175, 80, 0.3);
            display: inline-block;
            margin: 0 auto;
        }}
    </style>
    <div class="viewer-frame">
        {viewer_markup}
    </div>
    """
    components.html(framed_html, height=VIEW_HEIGHT + 80, width=VIEW_WIDTH + 80)


def set_show_flag():
    st.session_state["show_mol"] = True
    st.session_state.setdefault("queue_index", 0)


def move_queue(step: int) -> None:
    total = st.session_state.get("queue_len", 0)
    if not total:
        return
    st.session_state["queue_index"] = (st.session_state.get("queue_index", 0) + step) % total


def select_queue(idx: int) -> None:
    """Jump directly to a queued trajectory."""
    st.session_state["queue_index"] = idx
    st.session_state["show_mol"] = True


def extract_pdb_id(model_path: Path) -> str:
    """Derive the PDB id (e.g. 16pk_A) from a multi-model filename."""
    stem = model_path.stem
    suffix = "_multimodel"
    if stem.endswith(suffix):
        stem = stem[: -len(suffix)]
    parts = stem.split("_")
    if len(parts) >= 2:
        return "_".join(parts[:2])
    return stem




multimodel_files = discover_multimodel_files()
total_models = len(multimodel_files)
atlas_metadata = load_atlas_metadata()

if "queue_index" not in st.session_state:
    st.session_state["queue_index"] = 0

if total_models == 0:
    st.warning("No multi-model PDB files were found. Generate them before starting the viewer.")
else:
    st.session_state["queue_len"] = total_models
    st.session_state["queue_index"] = min(st.session_state["queue_index"], total_models - 1)

    set_show_flag()

    if st.session_state.get("show_mol", False):
        disable_navigation = total_models <= 1
        current_index = st.session_state.get("queue_index", 0)
        current_file = multimodel_files[current_index]
        current_label = current_file.relative_to(DATA_DIR).as_posix()
        pdb_id = extract_pdb_id(current_file)
        protein_info = atlas_metadata.get(pdb_id)

        with st.spinner("Loading multi-model trajectory..."):
            pdb_string = load_multimodel_pdb(str(current_file))

        viewer_col, info_col = st.columns([3, 2], gap="large")
        with viewer_col:
            render_mol(pdb_string)

        with info_col:
            if protein_info:
                st.write(f"**TITLE:** {protein_info['title'] or 'N/A'}")
                st.write(f"**UNIPROT ID:** {protein_info['uniprot'] or 'N/A'}")
                st.write(f"**ORGANISM:** {protein_info['organism'] or 'N/A'}")
                with st.expander("SEQUENCE", expanded=False):
                    st.text(protein_info["sequence"] or "No sequence found.")
            else:
                st.info("No ATLAS metadata found for this protein.")

            with st.container():
                st.subheader("AVAILABLE PROTEINS:")
                queue_labels: list[str] = []
                queue_index_map: dict[str, int] = {}
                for idx, model_path in enumerate(multimodel_files):
                    pdb_id = extract_pdb_id(model_path)
                    info = atlas_metadata.get(pdb_id, {})
                    label = info.get("uniprot") or pdb_id
                    title=info.get("title")
                    option_label = f"{idx + 1}. {label} ({title})"
                    queue_labels.append(option_label)
                    queue_index_map[option_label] = idx

                selected_label = st.radio(
                    "Select a trajectory",
                    queue_labels,
                    index=current_index,
                    key="queue_selector",
                )
                selected_idx = queue_index_map.get(selected_label, current_index)
                if selected_idx != current_index:
                    select_queue(selected_idx)
                    #st.experimental_rerun()

            nav_cols = st.columns(3)
            nav_cols[0].button(
                "Previous Model",
                on_click=move_queue,
                args=(-1,),
                disabled=disable_navigation,
            )
            nav_cols[2].button(
                "Next Model",
                on_click=move_queue,
                args=(1,),
                disabled=disable_navigation,
            )
