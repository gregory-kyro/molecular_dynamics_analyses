import MDAnalysis as mda
from MDAnalysis.analysis import align as mdaAlign
from MDAnalysis.analysis import rms

# define function to align trajectory
def align_traj(u, u_ref, segIDs, inMemory=True):
    u.trajectory[0]  # set first frame as reference for alignment
    alignment = mdaAlign.AlignTraj(u, u_ref, select="segid " + " ".join(segIDs) + " and not (name H* or name [123]H*)",
                                       verbose=True, in_memory=inMemory, weights="mass").run()
    return alignment
    
# create universe containing topology and trajectory
psf = 'path/to/protein/topology.pdb'
dcd = 'path/to/protein/trjaectory.dcd'
u = mda.Universe(psf,dcd, verbose=True, in_memory=True)
u_ref = mda.Universe('path/to/ref/topology.pdb','path/to/ref/trajectory.dcd', in_memory=True)
len(u.trajectory)

# align trajectory
align_traj(u, u_ref, segIDs='*', inMemory=True)

# compute RMSD
u.trajectory[-1]  # set mobile trajectory to last frame
u_ref.trajectory[0]  # set reference trajectory to first frame
mobile_ca = u.select_atoms('name CA')
ref_ca = u_ref.select_atoms('name CA')
rms.rmsd(mobile_ca.positions, ref_ca.positions, superposition=False)
