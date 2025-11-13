from pathlib import Path

import py3Dmol
import streamlit as st
import streamlit.components.v1 as components

DATA_DIR = Path("data")
VIEW_WIDTH = 500
VIEW_HEIGHT = 500

st.set_page_config(
    page_title="Phosphoglycerate Kinase Viewer",
    layout="wide",
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
            border: 2px solid #4caf50;
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


st.header("PHOSPHOGLYCERATE KINASE, GLYCOSOMAL")
st.write(
    "Queue up the generated multi-model trajectories below and animate them directly in the viewer."
)

multimodel_files = discover_multimodel_files()
total_models = len(multimodel_files)

if "queue_index" not in st.session_state:
    st.session_state["queue_index"] = 0

if total_models == 0:
    st.warning("No multi-model PDB files were found. Generate them before starting the viewer.")
else:
    st.session_state["queue_len"] = total_models
    st.session_state["queue_index"] = min(st.session_state["queue_index"], total_models - 1)

    st.button("Start Queue", on_click=set_show_flag)

    if st.session_state.get("show_mol", False):
        controls = st.columns(3)
        disable_navigation = total_models <= 1
        controls[0].button(
            "Previous Model",
            on_click=move_queue,
            args=(-1,),
            disabled=disable_navigation,
        )
        controls[2].button(
            "Next Model",
            on_click=move_queue,
            args=(1,),
            disabled=disable_navigation,
        )

        current_index = st.session_state.get("queue_index", 0)
        current_file = multimodel_files[current_index]
        current_label = current_file.relative_to(DATA_DIR).as_posix()

        st.subheader("Now Playing")
        st.write(f"**{current_label}** ({current_index + 1} of {total_models})")
        st.progress((current_index + 1) / total_models)

        with st.spinner("Loading multi-model trajectory..."):
            pdb_string = load_multimodel_pdb(str(current_file))
        render_mol(pdb_string)

        with st.expander("Queued trajectories", expanded=False):
            for idx, model_path in enumerate(multimodel_files):
                marker = ">>" if idx == current_index else "--"
                rel_path = model_path.relative_to(DATA_DIR).as_posix()
                st.write(f"{marker} {rel_path}")
