import os, pymol2

TOPO = "1j2l_A_analysis/1j2l_A.pdb"     # or .gro, .pdbqt, etc.
TRAJ = "1j2l_A_analysis/1j2l_A_R1.xtc"         # or .dcd/.trr
OUT_DIR = "frames"        # temporary frames
FPS = 3
QUALITY = 60              # 0-100

os.makedirs(OUT_DIR, exist_ok=True)

with pymol2.PyMOL() as pm:
    cmd = pm.cmd
    # speed up
    cmd.set("ray_shadows", 1)
    cmd.set("specular", 0)
    cmd.set("antialias", 0)
    cmd.set("cartoon_transparency", 1.0)


    # load topology and trajectory
    cmd.load(TOPO, "prot")
    cmd.load_traj(TRAJ, "prot", state=1, interval=1, average=0)

    # typical visualization
    cmd.hide("everything", "prot")
    cmd.show("cartoon", "prot")
    cmd.spectrum("count", "rainbow", "prot")
    cmd.orient("prot")

    # align frames to first frame for stable view (optional)
    cmd.intra_fit("prot")


    # after load_traj / intra_fit
    n_states = cmd.count_states("prot")
    assert n_states > 1, "Trajectory not loaded or only one state."

    # one-to-one mapping: state 1..N shown on frames 1..N
    cmd.mset(f"1 -{n_states}")       # NOT '1 x{n_states}'
    # optional: preview in GUI
    #cmd.mplay()
    #cmd.movie.roll(2, axis='z', first=1, last=180)
    #cmd.movie.produce("test_film.mp4", encoder="ffmpeg", quality=60)
    # render PNG frame sequence
    # use mpng to render all frames according to mset
    cmd.mpng(os.path.join(OUT_DIR, "frame"), width=512, height=328)