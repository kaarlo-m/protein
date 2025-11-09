import MDAnalysis as mda
import pyvista as pv

u = mda.Universe("data/1h02_B_analysis/1h02_B.pdb", "data/1h02_B_analysis/1h02_B_R1.xtc")
atoms = u.select_atoms("name CA")
sphere = pv.Sphere(radius=1.5)

# initial scene
points = pv.PolyData(atoms.positions)
glyphs = points.glyph(geom=sphere, scale=False)
plotter = pv.Plotter()
plotter.add_mesh(glyphs)

def update():
    for ts in u.trajectory:
        points.points = atoms.positions
        glyphs = points.glyph(geom=sphere, scale=False)
        plotter.clear()
        plotter.add_mesh(glyphs)
        plotter.render()

plotter.add_callback(update)
plotter.show()