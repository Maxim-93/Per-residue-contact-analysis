import MDAnalysis as mda
import random
import numpy as np
u = mda.Universe('r1_prot.gro','r1_prot.xtc')

#atom selection for c-alphas
all_coords=[]
for ts in u.trajectory:
    prot_CA=u.select_atoms('name CA')
    all_coords.append(prot_CA.positions)

transformed = []
# uncomment below line and re comment the next when test is complete
# for i in xrange(len(u.trajectory)):
# for i in xrange(len(prot_CA)):
for i in xrange(len(u.trajectory)):
    output_var=mda.lib.distances.distance_array(all_coords[i], all_coords[i])
    individual_density=np.where(output_var<=6, 1, 0)
    transformed.append(individual_density)

# establish the empty matrix that will contain the pre-normalised density matrix
density_matrix = []
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
    
print(density_matrix)
