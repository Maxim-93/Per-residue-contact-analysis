import collections
import numpy as np
from Bio import SeqIO
import os
import subprocess
import time
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

align_all = AlignIO.read("beta_a7_approved_MSA.clw", "clustal")
align_beta = AlignIO.read("beta_approved_MSA.clw", "clustal")
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
beta=subunit_orders(subunit="ACHB")

records= list(SeqIO.parse("beta_a7_approved_MSA.clw",'clustal'))
def reorder_subs(subunit_orders,subunit_name):
    os.remove('%s'%subunit_name+'.fasta')
    for i in subunit_orders:
        with open( '%s'%subunit_name+'.fasta',"a") as output_handle:
            SeqIO.write(records[i-1],output_handle,"fasta")
    return(list(SeqIO.parse('%s'%subunit_name+'.fasta','fasta')))

alpha7_block=reorder_subs(subunit_orders=alpha7,subunit_name="alpha7")
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
b_cons_dict=cons_dict(subunit=beta_block, order=beta)
# print(len(conservation_dictionary))

print(a7_cons_dict[2].most_common()[0][0])
print(len(conservation_dictionary))
def get_score(combined_alignment, subunit):
    test=[]
    for i in range(len(combined_alignment)):
        het_res=subunit[i].most_common()[0][0]
        het_freq=subunit[i].most_common()[0][1]
        comb_res=combined_alignment[i].most_common()[0][0]
        comb_freq=combined_alignment[i].most_common()[0][1]
        if het_res != comb_res and het_freq > comb_freq:
            test.append('not conserved')
        elif het_res != comb_res and het_freq < comb_freq:
            test.append('not conserved')
        elif het_res == comb_res:
            test.append('conserved')
    print(test)

get_score(combined_alignment=conservation_dictionary, subunit=b_cons_dict)
# # # This is a manual scoring function that you've made. Similar to kabat.
# def get_score(combined_alignment, subunit):
#     score=[]
#     combined_score=[]
#     for i in range(len(combined_alignment)):
#         # print((subunit[i].most_common()[0][1])/length)
#         score.append((subunit[i].most_common()[0][1]))
#         combined_score.append(conservation_dictionary[i].most_common()[0][1])
#     print(score)
#     print(combined_score)
#     # return(score)
#
# # alpha_score=get_score(combined_alignment=conservation_dictionary, subunit=a7_cons_dict, length=len(align_a7))
# beta_score=get_score(combined_alignment=conservation_dictionary, subunit=b_cons_dict)
