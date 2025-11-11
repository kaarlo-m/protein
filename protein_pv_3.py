#!/usr/bin/env python3
# cartoon_traj_pdbsec_pyvista.py
# Standalone PyVista viewer: plays an MD trajectory and renders secondary structure
# USING ONLY PDB HELIX/SHEET records (no DSSP). Helices=tubes, Sheets=ribbons+arrowheads, Coils=thin tubes.

import argparse, time, re
import numpy as np
import pyvista as pv
import MDAnalysis as mda
from MDAnalysis.analysis import align as mda_align

# ----------------------------- visual params -----------------------------
RAD_COIL    = 0.6     # tube radius for coils (Å)
RAD_HELIX   = 1.8     # tube radius for helices
SHEET_WIDTH = 4.0     # ribbon width for sheets
ARROW_LEN   = 4.0
ARROW_WIDTH = 4.0

CLR_HELIX = (0.85, 0.25, 0.25)   # RGB 0..1
CLR_SHEET = (0.25, 0.45, 0.85)
CLR_COIL  = (0.25, 0.75, 0.25)

# ----------------------------- PDB parsing -----------------------------
def parse_pdb_secondary(pdb_path):
    """
    Parse HELIX and SHEET records from a PDB file header.
    Returns two sets of residue keys: helices, sheets.
    Residue key: (chain_id, resseq:int, icode:str)
    """
    helices = set()
    sheets  = set()

    def key_range(chain_i, seq_i, ic_i, chain_e, seq_e, ic_e):
        # build inclusive residue range keys; insertion codes respected only at endpoints
        # For simplicity: include all integer seq between start and end; store with icode=''
        # and also include exact endpoints with given icode. Most trajectories drop icode anyway.
        keys = []
        if chain_i != chain_e:
            # Discontinuous across chains; treat as two separate single points
            seqs = [seq_i, seq_e]
            chains = [chain_i, chain_e]
            icodes = [ic_i, ic_e]
            for c, s, ic in zip(chains, seqs, icodes):
                keys.append((c, s, ic.strip() or ''))
                keys.append((c, s, ''))  # fallback key without icode
            return keys
        for s in range(min(seq_i, seq_e), max(seq_i, seq_e)+1):
            keys.append((chain_i, s, ''))
        # also add exact endpoints with insertion codes if present
        if ic_i.strip():
            keys.append((chain_i, seq_i, ic_i.strip()))
        if ic_e.strip():
            keys.append((chain_i, seq_e, ic_e.strip()))
        return keys

    with open(pdb_path, 'r', errors='ignore') as fh:
        for line in fh:
            rec = line[0:6]
            if rec.startswith('HELIX '):
                # PDB v3.3 format columns:
                # initChainID 20, initSeqNum 22-25, initICode 26
                # endChainID  32, endSeqNum  34-37, endICode  38
                try:
                    ci = line[19].strip() or ''
                    si = int(line[21:25])
                    ii = line[25].strip()
                    ce = line[31].strip() or ''
                    se = int(line[33:37])
                    ie = line[37].strip()
                    for k in key_range(ci, si, ii, ce, se, ie):
                        helices.add(k)
                except Exception:
                    continue
            elif rec.startswith('SHEET '):
                # initChainID 22, initSeqNum 23-26, initICode 27
                # endChainID  33, endSeqNum  34-37, endICode  38
                try:
                    ci = line[21].strip() or ''
                    si = int(line[22:26])
                    ii = line[26].strip()
                    ce = line[32].strip() or ''
                    se = int(line[33:37])
                    ie = line[37].strip()
                    for k in key_range(ci, si, ii, ce, se, ie):
                        sheets.add(k)
                except Exception:
                    continue
            elif rec.startswith('ATOM  ') or rec.startswith('HETATM'):
                # past header; stop scanning
                break
    return helices, sheets

def residue_key_from_mda_res(res):
    """
    Build keys to match PDB-derived keys for an MDAnalysis Residue.
    Returns a tuple: (primary_key, fallback_key)
      primary_key  = (chain_id, resid, icode) if available
      fallback_key = (chain_id, resid, '')    without insertion code
    Chain ID is taken from res.segid if present, otherwise from first atom's chainID if available; else ''.
    """
    # chain id
    chain = ''
    try:
        # MDAnalysis often stores PDB chain in .segid for PDB readers
        chain = (res.segid or '').strip()
    except Exception:
        chain = ''
    if not chain:
        try:
            chain_vals = res.atoms.chainIDs  # array if available
            if chain_vals is not None and len(chain_vals) > 0:
                chain = str(chain_vals[0]).strip()
        except Exception:
            chain = ''
    # resid
    resid = int(res.resid)
    # insertion code (not always available)
    icode = ''
    if hasattr(res, 'icode'):
        try:
            icode = (res.icode or '').strip()
        except Exception:
            icode = ''
    primary = (chain, resid, icode)
    fallback = (chain, resid, '')
    return primary, fallback

def build_static_ss_from_pdb(pdb_path, sel_residues):
    """
    Using HELIX/SHEET ranges from PDB header, return per-residue codes for selected residues:
    'H' helix, 'E' sheet, '-' coil.
    """
    helix_keys, sheet_keys = parse_pdb_secondary(pdb_path)
    ss_codes = []
    for res in sel_residues:
        primary, fallback = residue_key_from_mda_res(res)
        code = '-'
        if (primary in sheet_keys) or (fallback in sheet_keys):
            code = 'E'
        elif (primary in helix_keys) or (fallback in helix_keys):
            code = 'H'
        ss_codes.append(code)
    return np.array(ss_codes, dtype='<U1')

# ----------------------------- geometry builders -----------------------------
def contiguous_segments(ss_codes):
    if len(ss_codes) == 0:
        return
    cur = ss_codes[0]
    start = 0
    for i in range(1, len(ss_codes)):
        if ss_codes[i] != cur:
            yield (cur, np.arange(start, i))
            cur = ss_codes[i]
            start = i
    yield (cur, np.arange(start, len(ss_codes)))

def build_polyline(points_np):
    pd = pv.PolyData()
    pts = points_np.astype(np.float32, copy=False)
    pd.points = pts
    n = pts.shape[0]
    if n >= 2:
        cells = np.empty((n-1, 3), dtype=np.int32)
        cells[:, 0] = 2
        cells[:, 1] = np.arange(0, n-1, dtype=np.int32)
        cells[:, 2] = np.arange(1, n,   dtype=np.int32)
        pd.lines = cells
    return pd

def tube_from_polyline(polyline, radius):
    return polyline.tube(radius=radius, n_sides=16, capping=True)

def ribbon_from_polyline(polyline, width):
    return polyline.ribbon(width=width)

def arrow_head_at_end(polyline, length=ARROW_LEN, width=ARROW_WIDTH):
    pts = polyline.points
    if pts.shape[0] < 2:
        return None
    p_end  = pts[-1]
    p_prev = pts[-2]
    v = p_end - p_prev
    n = np.linalg.norm(v)
    if n == 0:
        return None
    dirv = v / n
    return pv.Cone(center=p_end - dirv * (length * 0.5),
                   direction=dirv, height=length, radius=width/2.0,
                   capping=True, resolution=24)

def compute_cartoon_meshes(ca_coords, ss_codes, add_arrows=True):
    out = []
    for typ, idxs in contiguous_segments(ss_codes):
        if len(idxs) == 0:
            continue
        if len(idxs) == 1:
            out.append((pv.Sphere(radius=RAD_COIL*1.2, center=ca_coords[idxs[0]]), CLR_COIL))
            continue
        poly = build_polyline(ca_coords[idxs])
        if typ == 'H':
            out.append((tube_from_polyline(poly, RAD_HELIX), CLR_HELIX))
        elif typ == 'E':
            rib = ribbon_from_polyline(poly, SHEET_WIDTH)
            out.append((rib, CLR_SHEET))
            if add_arrows:
                ah = arrow_head_at_end(poly)
                if ah is not None:
                    out.append((ah, CLR_SHEET))
        else:  # '-'
            out.append((tube_from_polyline(poly, RAD_COIL), CLR_COIL))
    return out

# ----------------------------- main -----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="PyVista trajectory viewer rendering helices/sheets from PDB HELIX/SHEET records (no DSSP).")
    ap.add_argument("-t", "--top", required=True, help="Topology PDB with HELIX/SHEET records")
    ap.add_argument("-x", "--traj", required=True, help="Trajectory: XTC/TRR/DCD/NC/…")
    ap.add_argument("--sel", default="name CA and protein", help="Selection (CA backbone recommended)")
    ap.add_argument("--align", action="store_true", help="Align all frames to first (CA)")
    ap.add_argument("--stride", type=int, default=1, help="Frame stride")
    ap.add_argument("--fps", type=int, default=24, help="Playback FPS")
    ap.add_argument("--window", type=int, nargs=2, default=[1280, 720], help="Window size W H")
    ap.add_argument("--no-arrows", action="store_true", help="Disable sheet arrowheads (faster)")
    ap.add_argument("--show-points", action="store_true", help="Also draw CA spheres/lines overlay")
    args = ap.parse_args()

    # Load trajectory
    u = mda.Universe(args.top, args.traj)
    sel_atoms = u.select_atoms(args.sel)
    if len(sel_atoms.residues) < 2:
        raise SystemExit("Selection too small; use at least CA backbone.")
    if u.filename != args.top:
        # MDAnalysis may read from in-memory; we still have the PDB path as args.top to parse header
        pass

    # Optional alignment
    if args.align:
        mda_align.AlignTraj(u, u, select="name CA and protein", in_memory=True).run()

    total_frames = len(u.trajectory)
    frames = list(range(0, total_frames, args.stride))
    if not frames:
        raise SystemExit("No frames after applying stride.")

    # Build static secondary structure codes per selected residue from PDB HELIX/SHEET
    ss_static = build_static_ss_from_pdb(args.top, sel_atoms.residues)  # array of 'H','E','-'

    # Initial frame
    u.trajectory[frames[0]]
    ca0 = sel_atoms.positions.copy()

    # Build initial meshes
    meshes0 = compute_cartoon_meshes(ca0, ss_static, add_arrows=not args.no_arrows)

    # Plotter and scene
    pv.set_plot_theme("default")
    pl = pv.Plotter(window_size=tuple(args.window))
    pl.set_background("white")
    pl.camera.zoom(1.25)

    # Optional overlays
    if args.show_points:
        pts_pd = pv.PolyData(ca0.astype(np.float32))
        pt_actor = pl.add_points(pts_pd, render_points_as_spheres=True, point_size=10.0, color=(0.1,0.1,0.1))
        ca_line = build_polyline(ca0)
        line_actor = pl.add_mesh(ca_line, color=(0.0,0.0,0.0), line_width=2.0)

    # Add cartoon actors
    actors = [pl.add_mesh(mesh, color=color, smooth_shading=True, opacity=1.0) for mesh, color in meshes0]

    # Playback state and controls
    state = {"i": 0, "paused": False, "running": True}

    def toggle_pause(): state.update(paused=not state["paused"])
    def step_fwd():     state.update(i=(state["i"] + 1) % len(frames)); update_frame()
    def step_back():    state.update(i=(state["i"] - 1) % len(frames)); update_frame()
    def quit_app():     state.update(running=False)

    pl.add_key_event("space", toggle_pause)
    pl.add_key_event("Right", step_fwd)
    pl.add_key_event("Left",  step_back)
    pl.add_key_event("q", quit_app)
    pl.add_key_event("Escape", quit_app)

    def rebuild_cartoon(ca_coords):
        nonlocal actors
        for a in actors:
            try: pl.remove_actor(a)
            except Exception: pass
        actors = []
        new_meshes = compute_cartoon_meshes(ca_coords, ss_static, add_arrows=not args.no_arrows)
        for mesh, color in new_meshes:
            actors.append(pl.add_mesh(mesh, color=color, smooth_shading=True, opacity=1.0))

    def update_frame():
        fi = frames[state["i"]]
        u.trajectory[fi]
        ca = sel_atoms.positions.copy()

        if args.show_points:
            pts_pd.points = ca.astype(np.float32)
            ca_line.points = ca.astype(np.float32)

        rebuild_cartoon(ca)
        pl.render()

    # Non-blocking show; drive loop manually
    pl.show(interactive_update=True, auto_close=False)

    target_dt = 1.0 / max(1, args.fps)
    last = time.perf_counter()

    while state["running"] and pl.ren_win is not None:
        pl.update()
        now = time.perf_counter()
        if not state["paused"] and (now - last) >= target_dt:
            update_frame()
            state["i"] = (state["i"] + 1) % len(frames)
            last = now
        time.sleep(0.001)

    pl.close()

if __name__ == "__main__":
    main()
