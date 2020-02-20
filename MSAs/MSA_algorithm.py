# Look for conserved sequences in obligate heteropentamers that are not conserved
# in obligate homopentamers. ID residues and validate against 3D structure.

# Look at sequence alignments of each subunit in each species.
# Find differences between conserved areas of heteromers where homomers are
# onserved. Then, compare conserved residues between subunit pairs to identify
# potential residues. Once you’ve ID’d potential residues- compare this
# (or “train” this) against 3D structure of two contacting subunits.
# Those residues that are conserved and within x distance will be deemed
# ‘important stoichiometric indicator’ (ISI).

from Bio.Seq import Seq
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
import collections

align_obligate_homomer = AlignIO.read("alpha7_approved_MSA.clw", "clustal")
align_gamma = AlignIO.read("gamma_approved_MSA.clw", "clustal")
align_delta = AlignIO.read("delta_approved_MSA.clw", "clustal")
align_epsilon = AlignIO.read("epsilon_approved_MSA.clw", "clustal")
align_all = AlignIO.read("approved_MSA.clw", "clustal")

# print(align_all[0])

# 1) load in alignments of both heteromer and obligate homomers
#  and combined alignments

# 2) get index ranges of each subunit and place them in individual groups

# As MSAs can become disordered- the following function will obtain the order your seqeunces occur in so that they can be
#  subsequently reordered.
def subunit_orders(subunit):
    orders=[]
    j=0
    for i in align_all:
        # print(i)
        align_string=str([i])
        # print(align_string)
        j=j+1
        if str(subunit) in align_string:
            orders.append(j)
    return(orders)

alpha7=subunit_orders(subunit="ACHA7")
delta=subunit_orders(subunit="ACHD")
epsilon=subunit_orders(subunit="ACHE")
gamma=subunit_orders(subunit="ACHG")

def reorder_subs(subunit_orders):
    subunit_block=[]
    for i in subunit_orders:
        subunit_block.append(align_all[i-1])
    return(subunit_block)

alpha7_block=reorder_subs(subunit_orders=alpha7)
delta_block=reorder_subs(subunit_orders=delta)
epsilon_block=reorder_subs(subunit_orders=epsilon)
gamma_block=reorder_subs(subunit_orders=gamma)

# slice through each column for a particular subunit block and calculate conservation
# for that point. Need to find a condition where if the character is an indel, give NA
# rather than an actual conservation score.
global_array=[]
for i in range(len(alpha7_block[0])):
    temp_array=[]
    for j in range(len(alpha7_block)):
        temp_array.append(align_all[j][i])
    # print(temp_array)
    global_array.append(temp_array)

for i in global_array:
    # print(i)
    frequencies = collections.Counter(i)
    print(frequencies)
# print(global_array)


# You need to update  this block of code to something more sophisticated.
# https://onlinelibrary.wiley.com/doi/full/10.1002/prot.10146#bib18
# contains information on conservation scores.
# What you want to do is make a function that uses the Zvelibil description.
# V_{zvelibil}=n_{const} * 1/10 (This is sort of a latex format for this equation)
# You need to take the 'truth table' and store it as a kind of Dictionary

# alpha7=align_all[0:homomer_length]
# delta=align_all[homomer_length:homomer_length+delta_length]
# epsilon=align_all[homomer_length+delta_length:homomer_length+delta_length+epsilon_length]
# gamma=align_all[homomer_length+homomer_length+delta_length:homomer_length+homomer_length+delta_length+epsilon_length]
# print(alpha7)
# print(delta)
# print(epsilon)
# print(len(gamma))
