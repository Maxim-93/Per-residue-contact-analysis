import MDAnalysis as mda
import random
import numpy as np
import matplotlib.pyplot as plt
import re

print('Enter beginning residue number of the second protein chain:')
x = input()


u = mda.Universe('del_beta.gro','del_beta_r1.xtc')

#atom selection for c-alphas
all_coords=[]
# Will store all CA atoms in a multidimensional array of trajectory length x number of CA atoms
for ts in u.trajectory:
    prot_CB=u.select_atoms('name CB')
    all_coords.append(prot_CB.positions)

# print(all_coords)

transformed = []
# Generates a new matrix for all contacts over every frame where distances are replaced with either a 0 or 1
# depending on if they meet the critera of being above or below 6 Angstroms.
for i in range(len(u.trajectory)):
    # MDAnalysis function to generate a per-residue distance matrix.
    output_var=mda.lib.distances.distance_array(all_coords[i], all_coords[i])
    # If the distance is below 8 angstroms then replace the distance with a 1.
    # If the distance is above 8 angstroms then replace the distance with a 0.
    individual_density=np.where(output_var<=8, 1, 0)
    # Append this transformed distance matrix to the array.
    transformed.append(individual_density)

print('transformed')

#
# # establish the empty matrix that will contain the pre-normalised density matrix
density_matrix = []
normalised_matrix = []
#This for loop is generating and i that corrosponds to the row to analyse
for i in range(len(transformed)):
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

print(len(density_matrix))

# # Normalise the density matrix
for i in range(len(prot_CB)):
    pre_norm_matrix=[]
    for j in range(len(prot_CB)):
        normalise=np.where(density_matrix[i][j]<=len(u.trajectory),round(float((float(density_matrix[i][j])/float(len(u.trajectory)))*100), 2),'blank')
        pre_norm_matrix.append(float(normalise))
    pre_norm_matrix=np.array(pre_norm_matrix)
    normalised_matrix.append(pre_norm_matrix)
normalised_matrix=np.array(normalised_matrix)


list_contacts=[]
prot=u.select_atoms('protein')
for i in range(0, int(x)-2):
    for j in range(int(x)-2, prot.n_residues-1):
        if normalised_matrix[i][j]>=70:
            list_contacts.append(i)
            list_contacts.append(j)

print("Delta residues")
delta_residues=[]
for i in list_contacts:
    if i < int(x)+1:
        delta_residues.append(i)
print(delta_residues)

print("Beta residues")
beta_residues=[]
for i in list_contacts:
    if i >= int(x):
        beta_residues.append(i)
print(beta_residues)

# # # Under construction:
# # # Give a table of residue pairs that are above a certain % cut-off.
# # # This will tell you what residues are important in forming contacts.
# # file=open("testfile.txt","w")
# # list_contacts=[]
# # for i in xrange(len(prot_CA)):
# #     sublist_contacts = []
# #     for j in xrange(len(prot_CA)):
# #         if normalised_matrix[i][j]>=95:
# #             file.write(str(i+1)+ " and " +str(j+1)+ "\n")
# # file.close()
#
# # Plotting information
# fig = plt.figure(figsize=(6, 3.2))
# ax = fig.add_subplot(111)
# ax.set_title('Contact Density')
# # plt.imshow(density_matrix)
# plt.imshow(normalised_matrix)
# ax.set_aspect('equal')
# cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
# cax.get_xaxis().set_visible(False)
# cax.get_yaxis().set_visible(False)
# cax.patch.set_alpha(0)
# cax.set_frame_on(False)
# # Set these axis to zoom in on different parts of the graph
# # ax.set_xlim(0, 250)
# # ax.set_ylim(250, 435)
# plt.colorbar(orientation='vertical')
# plt.show()
