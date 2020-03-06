import collections
import numpy as np
from Bio import SeqIO
import os
import subprocess
import time
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

# align_all = AlignIO.read("approved_MSA_2.clw", "clustal")
align_all = AlignIO.read("beta_a7_approved_MSA.clw", "clustal")
# align_gamma = AlignIO.read("gamma_approved_MSA.clw", "clustal")
align_beta = AlignIO.read("beta_approved_MSA.clw", "clustal")
# align_delta = AlignIO.read("delta_approved_MSA.clw", "clustal")
# align_epsilon = AlignIO.read("epsilon_approved_MSA.clw", "clustal")
align_a7=AlignIO.read("a7_approved_MSA.clw", "clustal")

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
# gamma=subunit_orders(subunit="ACHG")
# delta=subunit_orders(subunit="ACHD")
# epsilon=subunit_orders(subunit="ACHE")
beta=subunit_orders(subunit="ACHB")

records= list(SeqIO.parse("beta_a7_approved_MSA.clw",'clustal'))
def reorder_subs(subunit_orders,subunit_name):
    os.remove('%s'%subunit_name+'.fasta')
    for i in subunit_orders:
        with open( '%s'%subunit_name+'.fasta',"a") as output_handle:
            SeqIO.write(records[i-1],output_handle,"fasta")
    return(list(SeqIO.parse('%s'%subunit_name+'.fasta','fasta')))

alpha7_block=reorder_subs(subunit_orders=alpha7,subunit_name="alpha7")
# gamma_block=reorder_subs(subunit_orders=gamma,subunit_name="gamma")
# delta_block=reorder_subs(subunit_orders=delta,subunit_name="delta")
# epsilon_block=reorder_subs(subunit_orders=epsilon,subunit_name="epsilon")
beta_block=reorder_subs(subunit_orders=beta,subunit_name="beta")

global_array=[]
for i in range(len(align_all[0])):
    temp_array=[]
    for j in range(len(align_all)):
        temp_array.append(align_all[j][i])
    global_array.append(temp_array)
conservation_dictionary=[]
for i in global_array:
    frequencies = collections.Counter(i)
    conservation_dictionary.append(collections.Counter(i))

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
# g_cons_dict=cons_dict(subunit=gamma_block, order= gamma)
# d_cons_dict=cons_dict(subunit=delta_block, order=delta)
# e_cons_dict=cons_dict(subunit=epsilon_block, order=epsilon)
b_cons_dict=cons_dict(subunit=beta_block, order=beta)
print(len(a7_cons_dict))
print(len(b_cons_dict))
print(len(conservation_dictionary))
# This is the manual scoring function that you've made. Similar to kabat.
# for i in range(len(a7_cons_dict)):
#     print((a7_cons_dict[i].most_common()[0][1])/len(alpha7))

#     printa7
# print(a7_cons_dict)
# print(a7_cons_dict[i].most_common()[0][1])
# print(a7_cons_dict[2].most_common

# # os.remove('%s'%subunit+'.txt')
# conservation=[]
# # os.system("java -jar compbio-conservation-1.1.jar -i=%s_approved_MSA.clw -m=KABAT -o=%s.txt" % (subunit, subunit))
# os.system("java -jar compbio-conservation-1.1.jar -i=beta_a7_approved_MSA.clw -m=TAYLOR_GAPS -o=beta_a7.txt")

# with open("beta_a7.txt", "r") as file:
#     for line in file:
#         conservation.append(line)
# conservation=np.array(conservation)
# conservation=np.char.split(conservation)
# print(conservation)

def assign_score(subunit):
    os.remove('%s'%subunit+'.txt')
    conservation=[]
    # os.system("java -jar compbio-conservation-1.1.jar -i=%s_approved_MSA.clw -m=KABAT -o=%s.txt" % (subunit, subunit))
    os.system("java -jar compbio-conservation-1.1.jar -i=%s_approved_MSA.clw -m=ZVELIBIL -o=%s.txt" % (subunit, subunit))

    with open("%s.txt" % subunit, "r") as file:
        for line in file:
            conservation.append(line)
    conservation=np.array(conservation)
    conservation=np.char.split(conservation)
    return(conservation)

combined_conservation_score=assign_score(subunit='beta_a7')
a7_conservation =assign_score(subunit='alpha7')
# gam_conservation=assign_score(subunit='gamma')
# eps_conservation=assign_score(subunit='epsilon')
# del_conservation=assign_score(subunit='delta')
beta_conservation=assign_score(subunit='beta')

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
# gam_skip=skip_number(subunit_length=len(gamma_block), subunit=g_cons_dict)
# eps_skip=skip_number(subunit_length=len(epsilon_block), subunit=e_cons_dict)
# del_skip=skip_number(subunit_length=len(delta_block), subunit=d_cons_dict)
beta_skip=skip_number(subunit_length=len(beta_block), subunit=b_cons_dict)


p=1
for i in range(0,len(conservation_dictionary)):
    conservation_dictionary[i]['position']=i+1
    conservation_dictionary[i]['score']=combined_conservation_score[0][p]
    p=p+1


def local_append(dictionary, local_score, skip):
    p=1
    for i in range(0,len(dictionary)):
        dictionary[i]['position']=i+1
        if i not in skip:
            dictionary[i]['score']=local_score[0][p]
            p=p+1
    return(dictionary)

a7_local_dict= local_append(dictionary=a7_cons_dict,local_score=a7_conservation,  skip=a7_skip)
beta_local_dict=local_append(dictionary=b_cons_dict, local_score=beta_conservation, skip=beta_skip)
# print(a7_local_dict)
# print(beta_local_dict)
# beta_local_dict=local_append(dictionary=b_cons_dict, local_score=beta_conservation, skip=beta_skip)
# eps_local_dict=local_append(dictionary=e_cons_dict, local_score=eps_conservation, skip=eps_skip)
# gam_local_dict=local_append(dictionary=g_cons_dict, local_score=gam_conservation, skip=gam_skip)

# Give this function a rethink. The problem here is that if both homomer
# and heteromer are conserved at a specific point but not between them,
# then multiplication of both scores (1*1) will not successfully show differences between them.
# Simple multiplication of scores will only work here when one position is conserved
# and another is not. In order to get by this, obtain a global dictionary with both
# subunits. Count amino acids and append a score just like you did previously for individual
# subunits. You'll then want to make a comparison between the local score and this new global
# score.

# local_to_global_score=[]
# for i in range(len(conservation_dictionary)):
#     IF the heteromer is fully conserved (i.e 1) and
#     homomer is fully conserved (score 1), then the stoichiometry
#     indication score will also be 1.
#     Therefore, the lower the score the less likely
# print(conservation_dictionary[0]['score'])

# for i in range(len(conservation_dictionary):


print(conservation_dictionary)
print(beta_local_dict)
print(a7_local_dict)
# a7_local_dict position 1 score is not being shown for some reason...

# Calculate the difference between subunits by taking the difference between each subunits kabat conservation score ^2
# This transformation should help eccentuate differences between subunits where  large devaitions would indicate
# conserved-non-conserved pairs.

# The second simultaneous scoring funciton could then copae difference and square between the local and global alignments.
# This would tellyou if you have either conserved between subunit positions or conserved within subunits at that position.
def get_local_stoich_score(homomer, heteromer):
    score=[]
    compare_subs_dictionary={}
    for i in range(len(homomer)):
        if homomer[i]['score'] == 0 and heteromer[i]['score'] == 0:
            score.append('conserved gap')
            compare_subs_dictionary[i]='conserved gap'
        elif homomer[i]['score'] == 0 or heteromer[i]['score'] == 0:
            score.append('non conserved gap')
            compare_subs_dictionary[i]='non conserved gap'
        else:
            score.append((float(homomer[i]['score'])-float(heteromer[i]['score']))**2)
            # compare_subs_dictionary[i]['score']=str((float(homomer[i]['score'])-float(heteromer[i]['score']))**2)
            compare_subs_dictionary[i]=str((float(homomer[i]['score'])**2)-(float(heteromer[i]['score'])**2))
    # return(score)
    return(compare_subs_dictionary, score)

beta_stoich_score_dict, beta_stoich_score=get_local_stoich_score(homomer=a7_local_dict, heteromer=beta_local_dict)
print(beta_stoich_score_dict)

def get_global_stoich_score(global_alignment_score, heteromer_score):
    global_score_dict={}
    for i in range(0, len(global_alignment_score)):
        if heteromer_score[i]['score']=='non conserved gap':
            global_score_dict[i]='non conserved gap'
        elif heteromer_score[i]['score']=='conserved gap':
            global_score_dict[i]='non conserved gap'
        else:
            global_score_dict[i]=((float(global_alignment_score[i]['score'])**2)-(float(heteromer_score[i]['score'])**2))
            # print(heteromer_score[i])
    return(global_score_dict)

beta_global_score_dict=get_global_stoich_score(global_alignment_score=conservation_dictionary, heteromer_score=beta_local_dict)
print(beta_global_score_dict)
