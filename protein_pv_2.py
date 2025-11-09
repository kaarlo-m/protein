import argparse, time
import numpy as np
import pyvista as pv
import MDAnalysis as mda
from MDAnalysis.analysis import align as mda_align

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-t", "--top", required=True)
    ap.add_argument("-x", "--traj", required=True)
    ap.add_argument("--sel", default="name CA and protein")
    ap.add_argument("--align", action="store_true")
    ap.add_argument("--align-sel", default="name CA and protein")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--point-size", type=float, default=10.0)
    ap.add_argument("--window", type=int, nargs=2, default=[1280, 720])
    args = ap.parse_args()

    u = mda.Universe(args.top, args.traj)
    atoms = u.select_atoms(args.sel)
    n = len(atoms)

    line_cells = np.empty((n - 1, 3), dtype=np.int32)
    line_cells[:, 0] = 2              # number of points per segment
    line_cells[:, 1] = np.arange(0, n - 1)
    line_cells[:, 2] = np.arange(1, n)

    line_mesh = pv.PolyData()
    line_mesh.points = atoms.positions.astype(np.float32)
    line_mesh.lines = line_cells

    if args.align:
        mda_align.AlignTraj(u, u, select=args.align_sel, in_memory=True).run()

    # Build frame index honoring stride
    total = len(u.trajectory)
    frames = list(range(0, total, args.stride))
    if not frames:
        raise SystemExit("No frames (check inputs/stride).")

    # First frame
    u.trajectory[frames[0]]
    pts = pv.PolyData(atoms.positions.astype(np.float32))

    pl = pv.Plotter(window_size=tuple(args.window))
    pl.add_points(
        pts,
        render_points_as_spheres=True,
        point_size=args.point_size,
        opacity=1.0,
    )
    line_actor = pl.add_mesh(
        line_mesh,
        color="white",
        line_width=3,
        opacity=1.0,
    )
    pl.set_background("black")
    pl.camera.zoom(1.2)

    state = {"i": 0, "paused": False, "running": True}

    # Key bindings
    def toggle_pause():
        state["paused"] = not state["paused"]
    def step_forward():
        state["i"] = (state["i"] + 1) % len(frames)
        apply_frame()
    def step_back():
        state["i"] = (state["i"] - 1) % len(frames)
        apply_frame()
    def quit_app():
        state["running"] = False

    # Bind both "space" and literal space (backends differ)
    pl.add_key_event("space", toggle_pause)
    pl.add_key_event(" ", toggle_pause)
    pl.add_key_event("Right", step_forward)
    pl.add_key_event("Left", step_back)
    pl.add_key_event("q", quit_app)
    pl.add_key_event("Escape", quit_app)

    def apply_frame():
        u.trajectory[frames[state["i"]]]
        new_pos = atoms.positions.astype(np.float32)
        # Update spheres
        pts.points = new_pos
        # Update lines
        line_mesh.points = new_pos
        pl.render()

    # Open window without blocking; we’ll drive the loop
    pl.show(interactive_update=True, auto_close=False)

    target_dt = 1.0 / max(1, args.fps)
    last = time.perf_counter()

    while state["running"] and pl.ren_win is not None:
        now = time.perf_counter()
        # Process GUI events and keep window responsive
        pl.update()

        if not state["paused"] and (now - last) >= target_dt:
            apply_frame()
            state["i"] = (state["i"] + 1) % len(frames)
            last = now

        # Small sleep to avoid pegging a core
        time.sleep(0.001)

    pl.close()

if __name__ == "__main__":
    main()