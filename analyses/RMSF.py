import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSF

# create universe
a = mda.Universe(psf,dcd, verbose=True, in_memory=True)
a_ref = mda.Universe(psf_ref, dcd_ref, in_memory=True)
len(a.trajectory)

# select alpha carbon atoms
c1_alphas = a.select_atoms('protein and name CA')

# compute RMSF of alpha carbons
R1_rmsf= rms.RMSF(c1_alphas).run()

# add RMSF values as B-factors to pdb
a.add_TopologyAttr('tempfactors') # add empty attribute for all atoms
protein = a.select_atoms('protein') # select protein atoms
for residue, r_value in zip(protein.residues, R1_rmsf.rmsf):
    residue.atoms.tempfactors = r_value

# write universe as pdb file with RMSF as b factors for external visualization
aca = a.select_atoms('name CA')
with MDAnalysis.Writer('protein_RMSF.pdb', multiframe=True) as W:
    for ts in u.trajectory:
        W.write(aca)
