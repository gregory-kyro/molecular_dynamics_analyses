import MDAnalysis as mda
import MDAnalysis.analysis.pca as pca
from MDAnalysis.analysis import align as mdaAlign
import matplotlib.pyplot as plt

# create universe
prot_mobile=mda.Universe('top/path.pdb', 'traj/path.dcd', in_memory=True)

# align trajectory
align = mdaAlign

aligner1 = align.AlignTraj(prot_mobile, prot_mobile, select='name CA', in_memory=True).run()

# generate principal component matrices
pcu1 = pca.PCA(prot_mobile, select='protein and name CA',
             align=False, mean=None,
             n_components=None).run()

# print the number of alpha carbon atoms and the dimensions of the PC matrix
pc_1 = prot_mobile.select_atoms('protein and name CA')
n_pc1 = len(pc_1)
print('There are {} alpha carbon atoms in the WT analysis '.format(n_pc1))
print(pcu1.p_components.shape)

# print the variance of the first two principal components
print('The variance of the first principal component: {}'.format(pcu1.variance[0]))
print('The variance of the second principal component: {}'.format(pcu1.variance[1]))

# create plot of cumulative variance as a function of principal component
plt.plot(pcu1.cumulated_variance[:10])

# construct arrays of projection over first two principal components
transformed1 = pcu1.transform(pc_1, n_components=2)

# generate projected coordinates of first principal component 
pc1u1 = pcu1.p_components[:, 0]
trans1u1 = transformed1[:, 0]
projectedu1 = np.outer(trans1u1, pc1u1) + pcu1.mean
coordinatesu1 = projectedu1.reshape(len(trans1u1), -1, 3)

# create new universe
proj1 = mda.Merge(pc_1)
proj1.load_new(coordinatesu1, order="fac")

# write universe as xyz file for external visualization
proj1_ca = proj1.select_atoms('name CA')
with mda.Writer ('PC1.xyz',proj1_ca.n_atoms) as w:
  for ts in proj1.trajectory:
    w.write(proj1_ca)
    
# add pc rmsf as b factor
ca = proj1.select_atoms('protein and name CA')
pc1_rmsf = rms.RMSF(ca).run()
u = mda.Universe('/top/path.pdb','/traj/path.dcd')

u.add_TopologyAttr('tempfactors') # add empty attribute for all atoms
protein = u.select_atoms('protein') # select protein atoms
for residue, r_value in zip(protein.residues, pc1_rmsf.rmsf):
    residue.atoms.tempfactors = r_value

u_ca = u.select_atoms('name CA')
with MDAnalysis.Writer('PC1_RMSF.pdb', multiframe=True) as W:
    for ts in u.trajectory:
        W.write(u_ca)
