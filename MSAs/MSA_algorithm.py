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
import numpy as np
from Bio import SeqIO
import os
import subprocess
import time

align_obligate_homomer = AlignIO.read("alpha7_approved_MSA.clw", "clustal")
align_gamma = AlignIO.read("gamma_approved_MSA.clw", "clustal")
align_delta = AlignIO.read("delta_approved_MSA.clw", "clustal")
align_epsilon = AlignIO.read("epsilon_approved_MSA.clw", "clustal")
align_all = AlignIO.read("approved_MSA.clw", "clustal")

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

records= list(SeqIO.parse("approved_MSA.clw",'clustal'))
def reorder_subs(subunit_orders,subunit_name):
    os.remove('%s'%subunit_name+'.fasta')
    for i in subunit_orders:
        with open( '%s'%subunit_name+'.fasta',"a") as output_handle:
            SeqIO.write(records[i-1],output_handle,"fasta")
    return(list(SeqIO.parse('%s'%subunit_name+'.fasta','fasta')))

alpha7_block=reorder_subs(subunit_orders=alpha7,subunit_name="alpha7")
delta_block=reorder_subs(subunit_orders=delta,subunit_name="delta")
epsilon_block=reorder_subs(subunit_orders=epsilon,subunit_name="epsilon")
gamma_block=reorder_subs(subunit_orders=gamma,subunit_name="gamma")
# print(alpha7_block[1].seq)

# def conservation_score():
#     for index, record in enumerate(subunit_block):
#         print(record[0])
# Big problem here, when you try to run $java -jar compbio-conservation-1.1.jar -i=delta.fasta -m KABAT
#  in the command line, you get the error:
# "Input has badly aligned sequences with columns containing nothing but the gaps. Conservation methods cannot be calculated for such an alignment !"
# This means that you need to run this on a prealigned block, take the conservation score for that residue and map it back on to the corrosponding
# position that matches the position for the residue with respect to all the sequences. You might be able to do this with a python dictionary.

# Make a dictionary where you have a chronological list of residues with their corrosponding conservation score.

# slice through each column for a particular subunit block and calculate conservation
# for that point. Need to find a condition where if the character is an indel, give NA
# rather than an actual conservation score.
#
# Get conservation score, if Counter({'-': 6}) skip that iteration.
# if not Counter({'-': 6}), apply that conservation score to that dictionary element
global_array=[]
for i in range(len(alpha7_block[0])):
    temp_array=[]
    for j in range(len(alpha7_block)):
        temp_array.append(align_all[j][i])
    global_array.append(temp_array)

conservation_dictionary = []
for i in global_array:
    frequencies = collections.Counter(i)
    # print(frequencies)
    conservation_dictionary.append(collections.Counter(i))
print(conservation_dictionary)

def assign_score(subunit):
    # os.remove('%s'%subunit+'.txt')
    conservation=[]
    os.system("java -jar compbio-conservation-1.1.jar -i=%s_approved_MSA.clw -m=TAYLOR_GAPS -o=%s.txt" % (subunit, subunit))
    with open("%s.txt" % subunit, "r") as file:
        for line in file:
            conservation.append(line)
    conservation=np.array(conservation)
    conservation=np.char.split(conservation)
    return(conservation)

del_conservation=assign_score(subunit='delta')
a7_conservation =assign_score(subunit='alpha7')
eps_conservation=assign_score(subunit='epsilon')
gam_conservation=assign_score(subunit='gamma')

# determine what residues of the global alignment to skip over for that particular subbunit block
# when attributing a conservation score to the non-global alignment
def skip_number(subunit_length):
    gap_positions=[]
    h=0
    for i in conservation_dictionary:
        # print(h)
        elements=conservation_dictionary[h].elements()
        # print(elements)
        for value in elements:
            if value=='-' and i['-']==subunit_length:
                gap_positions.append(h)
        h=h+1
    gap_positions =list(dict.fromkeys(gap_positions))
    return(gap_positions)

a7_skip=skip_number(subunit_length=len(alpha7_block))
print(a7_skip)

# You need to update  this block of code to something more sophisticated.
# https://onlinelibrary.wiley.com/doi/full/10.1002/prot.10146#bib18
# contains information on conservation scores.
# What you want to do is make a function that uses the Zvelibil description.
# V_{zvelibil}=n_{const} * 1/10 (This is sort of a latex format for this equation)
# You need to take the 'truth table' and store it as a kind of Dictionary

# Actually want you can do is parse a list of FastaSequence objects as an input
# for the AACon command line method. You therefore need to find a way of saving
# your reordered subunit blocks as a FastaSequence object first.
