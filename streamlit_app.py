import streamlit as st
from stmol import showmol
import py3Dmol
import MDAnalysis as mda

pdb_path = "C:/Users/Kaarlo/protein/data/16pk_A_analysis/16pk_A.pdb"
xtc_path = "C:/Users/Kaarlo/protein/data/16pk_A_analysis/16pk_A_R1.xtc"  # not used directly by py3Dmol
u=mda.Universe(pdb_path, xtc_path)
ag=u.select_atoms("protein")
with mda.Writer("trajectory_multimodel.pdb", ag.n_atoms) as W:
    for ts in u.trajectory[::10]:  # stride frames to reduce file size
        W.write(ag)
pdb_str = open("trajectory_multimodel.pdb").read()

def render_mol(pdb_file):
    view = py3Dmol.view(width=500, height=500)
    view.addModelsAsFrames(pdb_str, "pdb")
    view.setStyle({"cartoon": {"color": "spectrum"}})
    view.animate({"interval": 75})
    view.setBackgroundColor('black')
    view.zoomTo()
    view.spin(True)
    #view.show()
    showmol(view, height=500, width=500)
    

def set_show_flag():
    st.session_state["show_mol"] = True



_ = st.sidebar.text_area('Show the protein')
st.sidebar.button('Show', on_click=set_show_flag)
st.header("PHOSPHOGLYCERATE KINASE, GLYCOSOMAL")
if st.session_state.get("show_mol", False):
    render_mol(pdb_path)