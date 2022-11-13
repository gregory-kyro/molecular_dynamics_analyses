# graph_theoretic_protein_analyses -- Conceptual Overview

## ROOT MEAN SQUARE DEVIATION (RMSD) - OVERVIEW

RMSD measures the deviation of a target set of coordinates (i.e. a structure) to a reference set of coordinates, with RMSD=0 indicating a perfect overlap. Then it follows that if we have a MD trajectory one would expect that the lower the RMSD the less changes happen in the time scale studied. 

RMSD is defined as:

<p align="center">
$RMSD = \sqrt{\frac{\sum_{i} m_i(X_i - Y_i)^2}{M}}$
</p>

where **N** is the number of atoms,  $m_{i}$ is the mass of atom $i$, $X_i$ is the coordinate vector for target atom $i$, $Y_i$ is the coordinate vector for reference atom $i$, and $M$ is the total mass. If the $RMSD$ is not mass-weighted, for all $i$, $m_i = 1$, and $M = N$.
When calculating RMSD of a target to reference structure, there are two very important requirements:
1. The number of atoms in the target must match the number of atoms in the reference.
2. The ordering of atoms in the target must match the ordering of atoms in the reference.

## ROOT MEAN SQUARE FLUCTUATION (RMSF) - OVERVIEW

The root-mean-square-fluctuation (RMSF) of a structure is the time average of the RMSD. It is calculated according to the below equation, where $x_i$ is the coordinates of atom $i$ and $⟨x_i⟩$ is the ensemble average position of $i$:

<p align="center">
$RMSF_{(i)}=\sqrt{\frac{1}{N}\sum_{i}{⟨(x_i−⟨x_i⟩)^2⟩}}$
</p>

Where the RMSD quantifies how much a structure diverges from a reference over time, the RMSF can reveal which areas of the system are the most mobile. While RMSD is frequently calculated to an initial state, the RMSF should be calculated to an average structure of the simulation. An area of the structure with high RMSF values frequently diverges from the average, indicating high mobility. When RMSF analysis is carried out on proteins, it is typically restricted to backbone or alpha-carbon atoms; these are more characteristic of conformational changes than the more flexible side-chains.

## PRINCIPAL COMPONENT ANALYSIS (PCA) - OVERVIEW

PCA decomposes a system of observations into linearly uncorrelated variables called principal components. The components are ordered; the first component accounts for the largest variance in the data, with following components accounting for progressively lower variance. PCA is often applied to MD trajectories to extract large-scale conformational motions (i.e. 'essential dynamics') of a protein. Greater variance = greater conformational motions. The frame-by-frame conformational fluctuation can be considered a linear combination of the essential dynamics yielded by the PCA.

Steps:
 1. Align each frame in trajectory.
 2. Contruct a $3N * 3N$ covariance matrix for the *N* atoms in trajectory.     
 3. Diagonalize covariance matrix. Principal components are the eigenvectors 
    of the covariance matrix. The associated variance are the eigenvalues.
 4. Sort the eigenvalues (variance) so that the principal components 
    (eigenvectors) are ordered by variance. 
