import MDAnalysis as mda
from MDAnalysis.coordinates.XTC import XTCWriter
from MDAnalysis.analysis import align
from MDAnalysis.transformations import unwrap, center_in_box, wrap
import sys

pdb_in = "../4_openmm/step4_equilibration.pdb"   # PDB used when running OpenMM (must match DCD atom order)
dcds = ["../mhh1/step5_1.dcd", "../mhh1/step5_2.dcd"]
        
pdb_out = "mhhr1_clean.pdb"
xtc_out = "mhhr1_clean_fit.xtc"
reference = pdb_in      # ref for fitting (usually the initial PDB)

# Load Universe 
u = mda.Universe(pdb_in, dcds)

# Select protein only, no drude & solvent
prot = u.select_atoms("protein and not (name D* or name DC* or name LP*)")

print(f"Original atom count: {len(u.atoms)}, kept atoms: {len(prot)}")

# Add transformations: unwrap, center, wrap
u.trajectory.add_transformations(
    unwrap(prot),
    center_in_box(prot, center='geometry'),
    wrap(prot, compound='residues'),
)

# Write cleaned PDB
prot.write(pdb_out)
print(f"Wrote cleaned PDB: {pdb_out}")

# Align traj on protein backbone
ref = mda.Universe(reference)   # reference for alignment
aligner = align. AlignTraj(u, ref, 
			select="backbone and not (name D* or name DC* or name LP)", 
			in_memory=True) 
aligner.run() 

# Write aligned traj
with XTCWriter(xtc_out, prot.n_atoms) as W:
    for ts in u.trajectory:
        W.write(prot)
print(f"Wrote cleaned fitted XTC: {xtc_out}")


