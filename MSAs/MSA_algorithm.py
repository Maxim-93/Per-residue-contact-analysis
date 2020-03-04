
# Look for conserved sequences in obligate heteropentamers that are not conserved
# in obligate homopentamers. ID residues and validate against 3D structure.
# Look at sequence alignments of each subunit in each species.
# Find differences between conserved areas of heteromers where homomers are
# onserved. Then, compare conserved residues between subunit pairs to identify
# potential residues. Once you’ve id’d potential residues- compare this
# (or “train” this) against 3D structure of two contacting subunits.
# Those residues that are conserved and within x distance will be deemed
# ‘important stoichiometric indicator’ (ISI).
import collections
import numpy as np
from Bio import SeqIO
import os
import subprocess
import time
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

align_obligate_homomer = AlignIO.read("alpha7_approved_MSA.clw", "clustal")
align_delta = AlignIO.read("delta_approved_MSA.clw", "clustal")
align_epsilon = AlignIO.read("epsilon_approved_MSA.clw", "clustal")
align_beta = AlignIO.read("beta_approved_MSA.clw", "clustal")
align_gamma = AlignIO.read("gamma_approved_MSA.clw", "clustal")
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
        align_string=str([i])
        j=j+1
        if str(subunit) in align_string:
            orders.append(j)
    return(orders)

alpha7=subunit_orders(subunit="ACHA7")
delta=subunit_orders(subunit="ACHD")
epsilon=subunit_orders(subunit="ACHE")
gamma=subunit_orders(subunit="ACHG")
beta=subunit_orders(subunit="ACHB")

records= list(SeqIO.parse("approved_MSA.clw",'clustal'))
def reorder_subs(subunit_orders,subunit_name):
    # os.remove('%s'%subunit_name+'.fasta')
    for i in subunit_orders:
        with open( '%s'%subunit_name+'.fasta',"a") as output_handle:
            SeqIO.write(records[i-1],output_handle,"fasta")
    return(list(SeqIO.parse('%s'%subunit_name+'.fasta','fasta')))

alpha7_block=reorder_subs(subunit_orders=alpha7,subunit_name="alpha7")
delta_block=reorder_subs(subunit_orders=delta,subunit_name="delta")
epsilon_block=reorder_subs(subunit_orders=epsilon,subunit_name="epsilon")
gamma_block=reorder_subs(subunit_orders=gamma,subunit_name="gamma")
beta_block=reorder_subs(subunit_orders=beta,subunit_name="beta")

print(beta_block)

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

def cons_dict(subunit, order):
    global_array=[]
    #for each column, go through every row for a particular subunit and append that row to an array which can later be sorted by
    #collections.Counter
    for i in range(len(subunit[0])):
        temp_array=[]
        for j in order:
            temp_array.append(align_all[j-1][i])
        global_array.append(temp_array)
    conservation_dictionary = []
    for i in global_array:
        frequencies = collections.Counter(i)
        conservation_dictionary.append(collections.Counter(i))
    return(conservation_dictionary)

a7_cons_dict=cons_dict(subunit=alpha7_block, order=alpha7)
d_cons_dict=cons_dict(subunit=delta_block, order=delta)
e_cons_dict=cons_dict(subunit=epsilon_block, order=epsilon)
g_cons_dict=cons_dict(subunit=gamma_block, order=gamma)
b_cons_dict=cons_dict(subunit=beta_block, order=beta)

# This block utilises the Barton group AACon program.
# Try to put your own kabat function in rather than relying on AACon in the future.

# VKabat= k/n1 * N
# where k is the number of amino types present at the
# aligned position, n1 is the number of times the most
# commonly occurring amino acid appears there, and N is
# the number of sequences in the alignment. The variable N
# acts as a scaling factor and is constant for a given
# alignment. (PROTEINS: Structure, Function, and Genetics 48:227–241 (2002))

def assign_score(subunit):
    os.remove('%s'%subunit+'.txt')
    conservation=[]
    os.system("java -jar compbio-conservation-1.1.jar -i=%s_approved_MSA.clw -m=KABAT -o=%s.txt" % (subunit, subunit))

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
beta_conservation=assign_score(subunit='beta')

# determine what residues of the global alignment to skip over for that particular subunit block
# when attributing a conservation score to the non-global alignment.
def skip_number(subunit_length, subunit):
    gap_positions=[]
    h=0
    for i in subunit:
        elements=subunit[h].elements()
        for value in elements:
            if value=='-' and i['-']==subunit_length:
                gap_positions.append(h)
        h=h+1
    gap_positions =list(dict.fromkeys(gap_positions))
    return(gap_positions)

a7_skip=skip_number(subunit_length=len(alpha7_block), subunit=a7_cons_dict)
eps_skip=skip_number(subunit_length=len(epsilon_block), subunit=e_cons_dict)
del_skip=skip_number(subunit_length=len(delta_block), subunit=d_cons_dict)
gam_skip=skip_number(subunit_length=len(gamma_block), subunit=g_cons_dict)
# beta_skip=skip_number(subunit_length=len(beta_block), subunit=b_cons_dict)

# As ?_conservation dictionary (contains the scores) only considers local alignment,
# you need to go through the ?_conservation dictionaries and append that score to the
# ?_cons_dict, adding a new dictionary variable 'score' and skipping this appending
# if ?_skip number comes up. Note that ?_conservation
# dictionary index actually starts at 1 as the 0th index is just 'TAYLOR_GAPS'
# or whatever score you use in the AACon step.
# You also need to append the global conservation score to this dictionary as well.

def local_append(dictionary, local_score, skip):
    p=1
    for i in range(0,len(dictionary)):
        if i not in skip:
            dictionary[i]['local score']=local_score[0][p]
            p=p+1
    return(dictionary)

a7_local_dict= local_append(dictionary=a7_cons_dict,local_score=a7_conservation,  skip=a7_skip )
del_local_dict=local_append(dictionary=d_cons_dict, local_score=del_conservation, skip=del_skip)
# beta_local_dict=local_append(dictionary=b_cons_dict, local_score=beta_conservation, skip=beta_skip)
eps_local_dict=local_append(dictionary=e_cons_dict, local_score=eps_conservation, skip=eps_skip)
# gam_local_dict=local_append(dictionary=g_cons_dict, local_score=gam_conservation, skip=gam_skip)

# print(del_local_dict)

# Gamma doesn't seem to work for some reason for local_append function.
# Now that you have the local scores for each subunit, you need to
# compare the local score of a heteromer with a homomer and output
# another score, indicating the potential stoichiometric importance
# of that residue.

# As all gaps wont have a score, you won't be able to calculate an output
# Need to figure a way around this. See followign lines for example.
# [1,2,7,4,-,1,2,9,5]
# [2,-,2,1,4,7,8,5,-]

# I think this function is slightly flawed. Have a rethink....
def stoich_score(homomer, heteromer):
    score=[]
    for i in range(len(homomer)):
        if homomer[i]['local score'] == 0 and heteromer[i]['local score'] == 0:
            score.append('conserved gap')
        elif homomer[i]['local score'] == 0 or heteromer[i]['local score'] == 0:
            score.append('non conserved gap')
        else:
            score.append(float(homomer[i]['local score'])*float(heteromer[i]['local score']))
    return(score)

del_stoich_score=stoich_score(homomer=a7_local_dict, heteromer=del_local_dict)
eps_stoich_score=stoich_score(homomer=a7_local_dict, heteromer=eps_local_dict)

# beta_stoich_score=stoich_score(homomer=a7_local_dict, heteromer=beta_local_dict)
# print(del_stoich_score)
# print(eps_stoich_score)
# print(beta_stoich_score)

# You need to update  this block of code to something more sophisticated.
# https://onlinelibrary.wiley.com/doi/full/10.1002/prot.10146#bib18
# contains information on conservation scores.
# What you want to do is make a function that uses the Zvelibil description.
# V_{zvelibil}=n_{const} * 1/10 (This is sort of a latex format for this equation)
# You need to take the 'truth table' and store it as a kind of Dictionary

# Actually want you can do is parse a list of FastaSequence objects as an input
# for the AACon command line method. You therefore need to find a way of saving
# your reordered subunit blocks as a FastaSequence object first.
