# This script will be more detailed than our other contact analysis matrix insofar that instead of just picking CA atoms,
# For each timestep, go through each individual residue of the protein and identify what amino acid it is. Then select the residues
# R group and use that to calculate distances. If these residues occur at interfaces, then it is likely that it's the R groups that
# will be contributing to these interface contacts.

# Another analysis to do might be to identify individual protein chains and, if a residue's R group falls within X distance of another selected
# protein chain, get the resnumber.

import MDAnalysis as mda
import random
import numpy as np
import matplotlib.pyplot as plt
import re
u = mda.Universe('test.gro','r1_prot.xtc')

# A word of caution... Make sure that your .gro file protein numbering is
# set so that resnumbers don't skip.
# If in doubt, $gmx editconf -f '?.gro' -resnr 1 -o whatever.gro

#Dictionary with corrosponding R groups for residue names
amino_acid_r_names = {
    'ARG':'CZ',
    'HIS':'CE1',
    'LYS':'NZ',
    'ASP':'CG',
    'GLU':'CD',
    'SER':'OG',
    'THR':'OG1',
    'ASN':'CG',
    'GLN':'CD',
    'CYS':'SG',
    'GLY':'CA',
    'PRO':'CG',
    'ALA':'CB',
    'VAL':'CB',
    'ILE':'CD',
    'LEU':'CG',
    'MET':'SD',
    'PHE':'CZ',
    'TYR':'OH',
    'TRP':'NE1'
}

all_coords=[]
for ts in u.trajectory[0:2]:
# for ts in u.trajectory:
    prot=u.select_atoms('protein')
    frame_specific_array=[]
    for j in range(1, prot.n_residues):
        # Go through each residue chronologically according to the for loop and get
        # that residues resname
        residue=prot.select_atoms('resid %i' % j).resnames[0]
        #Access the name from the Dictionary
        R_name=amino_acid_r_names.get(residue)
        residue_to_append=u.select_atoms('resid %i and name %s' % (j, R_name))
        position=residue_to_append.positions
        # Extend rather than append in order to merge the lists so that we have something like this:
        # list=['var1', 'var2', 'var3', 'var4', 'var5'] rather than
        # list=['var1', 'var2', 'var3', ['var4', 'var5']]
        frame_specific_array.extend(position)
    all_coords.append(frame_specific_array)

# all_coords is an array ordered with each element containing the positions of all the residues for each timestep.
# This gives it dimensions of max(prot.n_residues) x number of trajectory frames specified
all_coords=np.array(all_coords)
print(len(all_coords[0]))
# print(all_coords)

transformed = []
# Generates a new matrix for all contacts over every frame where distances are replaced with either a 0 or 1
# depending on if they meet the critera of being above or below 6 Angstroms.
# Go through the distances for each frame in trajectory. Note that each frame of the trajectory should be equivalent to len(all_coords)
for i in range(0,len(all_coords)):
# for i in range(len(u.trajectory)):
    # MDAnalysis function to generate a per-residue distance matrix.
    # This is a matrix of dimensions max(prot.n_residues) x max(prot.n_residues)
    output_var=mda.lib.distances.distance_array(all_coords[i], all_coords[i])
    # If the distance is below 6 angstroms then replace the distance with a 1.
    # If the distance is above 6 angstroms then replace the distance with a 0.
    individual_density=np.where(output_var<=6, 1, 0)
    # print(individual_density)
    # Append this transformed distance matrix to the array.
    transformed.append(individual_density)
    print(len(output_var[0]))
    # print(output_var)

# print(len(transformed[0]))
# establish the empty matrix that will contain the pre-normalised density matrix
density_matrix = []
normalised_matrix = []
#This for loop is generating and i that corrosponds to the row to analyse
# for i in range(1, prot.n_residues):
for i in range(0, len(transformed[0])):
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

# print(len(transformed[0]))
# print(density_matrix)
# Normalise the density matrix
for i in range(0, prot.n_residues-1):
    pre_norm_matrix=[]
    for j in range(0, prot.n_residues-1):
        normalise=np.where(density_matrix[i][j]<=len(transformed),round(float((float(density_matrix[i][j])/float(len(transformed)))*100), 2),'blank')
        pre_norm_matrix.append(float(normalise))
    pre_norm_matrix=np.array(pre_norm_matrix)
    normalised_matrix.append(pre_norm_matrix)
normalised_matrix=np.array(normalised_matrix)


# # # Under construction:
# # # Give a table of residue pairs that are above a certain % cut-off.
# # # This will tell you what residues are important in forming contacts.


# Go through normalised matrix, if normalised_matrix[i][j] => 70 , append i+1 and j+1 to a list
# Ok, this works but there is alot of redundancy in the list. Must figure out how to remove repeated and self-matching pairs of residues.
# Also, perhaps rather than listing the pairs explicitly- just list the amino acids that make these contacts themselves.
# file=open("testfile.txt","w")
# list_contacts=[]
# for i in range(0, prot.n_residues-1):
#     sublist_contacts = []
#     for j in range(0, prot.n_residues-1):
#         if normalised_matrix[i][j]>=70:
#             file.write(str(i+1)+ " and " +str(j+1)+ "\n")
#
# file.close()

# For this to work properly, you need to specify what reisdues corrospond to each chain. This way, you can remove all residues that interact with
# other residues in the same chain.
list_contacts=[]

# Create an input for where chain 2 starts.
print('Enter beginning residue number of the second protein chain:')
x = input()
# For test.gro, this is 213
for i in range(0, int(x)-2):
    for j in range(int(x)-2, prot.n_residues-1):
        if normalised_matrix[i][j]>=70:
            list_contacts.append(i)
            list_contacts.append(j)

print(set(list_contacts))

# As MacPyMOL embeds Python directly, you cannot take this list of contacts and put them into pymol within
# this script. You will need to print them out and place them into a .pml file.

# Plotting information
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
