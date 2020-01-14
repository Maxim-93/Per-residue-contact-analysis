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

density_matrix = []
# for i in xrange(len(prot_CA)):
for i in xrange(len(prot_CA)):
    calc_matrix=[]
    for frame in transformed:
        calc_matrix.append(frame[i])
        # print(len(frame))
        # print(frame)
    # print(calc_matrix)
    calc_matrix=np.array(calc_matrix)
    row=calc_matrix.sum(axis=0)
    density_matrix.append(row)
    
print(density_matrix)
