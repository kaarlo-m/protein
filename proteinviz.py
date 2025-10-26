import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD
import nglview as nv

# Load simulation results with a single line
u = mda.Universe(PSF, DCD)

# Select atoms
ag = u.select_atoms('name OH')

# Atom data made available as Numpy arrays
ag.positions
#ag.velocities
#ag.forces

# Iterate through trajectories
for ts in u.trajectory:
    print(ag.center_of_mass())

view=nv.show_mdanalysis(ag)
view