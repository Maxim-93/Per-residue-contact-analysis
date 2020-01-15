import MDAnalysis as mda
import random
import numpy as np
import matplotlib.pyplot as plt
import re
u = mda.Universe('r1_prot.gro','r1_prot.xtc')

#atom selection for c-alphas
all_coords=[]
for ts in u.trajectory:
    prot_CA=u.select_atoms('name CA')
    all_coords.append(prot_CA.positions)

transformed = []
# Generates a new matrix for all contacts over every frame where distances are replaced with either a 0 or 1
# depending on if they meet the critera of being above or below 6 Angstroms.
for i in xrange(len(u.trajectory)):
    # MDAnalysis function to generate a per-residue distance matrix.
    output_var=mda.lib.distances.distance_array(all_coords[i], all_coords[i])
    # If the distance is below 6 angstroms then replace the distance with a 1.
    # If the distance is above 6 angstroms then replace the distance with a 0.
    individual_density=np.where(output_var<=10, 1, 0)
    # Append this transformed distance matrix to the array.
    transformed.append(individual_density)

# establish the empty matrix that will contain the pre-normalised density matrix
density_matrix = []
normalised_matrix = []
#This for loop is generating and i that corrosponds to the row to analyse
for i in xrange(len(prot_CA)):
    # This is a temporary matrix where you store the i-th row for every consecutive frame. Each element in these rows will then
    # be added to make one row, before then being refreshed for storing the next row in the series.
    calc_matrix=[]
    for frame in transformed:
        # Append the i-th row of all frames to 'calc_matrix'
        calc_matrix.append(frame[i])
    # make 'calc_matrix' a numpy array so that it can be summed
    calc_matrix=np.array(calc_matrix)
    # add each connected element in every row and save to make a single 'row'
    row=calc_matrix.sum(axis=0)
    # append this variable to the density matrix
    density_matrix.append(row)

# Normalise the density matrix
for i in xrange(len(prot_CA)):
    pre_norm_matrix=[]
    for j in xrange(len(prot_CA)):
        normalise=np.where(density_matrix[i][j]<=len(u.trajectory),round(float((float(density_matrix[i][j])/float(len(u.trajectory)))*100), 2),'blank')
        pre_norm_matrix.append(float(normalise))
    pre_norm_matrix=np.array(pre_norm_matrix)
    normalised_matrix.append(pre_norm_matrix)
normalised_matrix=np.array(normalised_matrix)

# Plotting information
fig = plt.figure(figsize=(6, 3.2))
ax = fig.add_subplot(111)
ax.set_title('Contact Density')
# plt.imshow(density_matrix)
plt.imshow(normalised_matrix)
ax.set_aspect('equal')
cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)
# Set these axis to zoom in on different parts of the graph
# ax.set_xlim(0, 250)
# ax.set_ylim(250, 435)
plt.colorbar(orientation='vertical')
plt.show()
