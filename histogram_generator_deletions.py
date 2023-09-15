#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from Bio import SeqIO
import argparse
import multiprocessing
import re
from skbio.alignment import StripedSmithWaterman
import os
import csv
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import sys

#Parse input
#_____________________________________________________________________________________________________________________________
parser = argparse.ArgumentParser(description=""" --- """)

parser.add_argument('-T', '--threads', type=int, default=0, help='Number of threads to execute in parallel. Defaults to CPU count.')
parser.add_argument('--variant_merge_file', help='Output "variant merge file" .tsv from dimsum pipeline', required=True)
parser.add_argument('--round', help='Enter selection round: either "r1" or "r2"', required=True)

args = parser.parse_args()
threads = args.threads
if threads == 0:
    threads = multiprocessing.cpu_count()
    
#_____________________________________________________________________________________________________________________________
#set up input variables
variant_merge_file = args.variant_merge_file
output_round = str(args.round)

output_round_options = ["r1", "r2"]
if output_round not in output_round_options:
    print("round not correctly specified: enter einter r1 or r2")
    sys.exit()

variant_merge_list = []
#fills up dictionary
with open (variant_merge_file, 'rt') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader: 
        variant_merge_list.append(row)
#get rid of first row for analysis
variant_merge_list.pop(0)

#get sequence counts (ci) from raw data list
sequences = [str(x[0]) for x in variant_merge_list]
c_input_list = [int(x[8]) for x in variant_merge_list]
c1_list = [int(x[9]) for x in variant_merge_list]
c2_list = [int(x[10]) for x in variant_merge_list]

#_____________________________________________________________________________________________________________________________
Nt_input = 0
Nt_1 = 0
Nt_2 = 0

for i in range(0, len(c_input_list)):
    Nt_input += c_input_list[i]
    Nt_1 += c1_list[i]
    Nt_2 += c2_list[i]

#define total number based on round
if output_round == "r1":
    Nt = Nt_1
if output_round == "r2":
    Nt = Nt_2

#reference:
reference = "NNNCATAAAAACVRCVRCRRCRRCAAGGATAACNDTCATGATNDTGATAACCATCTGCAGAACGTGATTGAAGATATTCATGATTTTATGCAGGGCGGCGGCAGCGGCGGCAAACTGCAGGAAATGATGAAAGAATTCCAGCAGGTGCTGGATGAANDTAACAACVRCVRCVRCRRCRRCAAACATNDTNDT"

#_____________________________________________________________________________________________________________________________
#alignment function

def aln_target(target, reference):
    alignments = []
    target_sequences = [target]
    query = StripedSmithWaterman(reference)
    
    for target_sequence in target_sequences:
        alignment = query(target_sequence)
        alignments.append(alignment)
    target_aln = alignments[0].aligned_target_sequence
    return(target_aln)

#_____________________________________________________________________________________________________________________________
#generate data frames and csv files for histograms

data_frame_positions = list(np.arange(4, 193, 1))
data_frame_initial_fill = []
for i in range(0, 193):
    data_frame_initial_fill.append(0)

df_input_del_positions_counts = pd.DataFrame(list(zip(data_frame_positions, data_frame_initial_fill)), columns = ['position', 'count'])
df_output_del_positions_counts = pd.DataFrame(list(zip(data_frame_positions, data_frame_initial_fill)), columns = ['position', 'count'])

for i in range(0, len(sequences)):
    sequence = sequences[i]
    alignment = aln_target(sequence, reference)
    #exclude more that one deletion
    count_del = 0
    for a in range(len(alignment)):
        if alignment[a] == "-":
            count_del += 1
    #if alignment starts with constant region -> discards truncated sequences, and sequences with more than one deletion
    if alignment[0:3] == "cat" and count_del < 2:
        for n in range(len(alignment)):
            #loop through alignment, excluding more than single nucleotide deletions:
            if alignment[n] == "-" and n < 193:
                if output_round == "r1":
                    #get old value in dataframe, add current count and then replace the old with the new one
                    old = df_output_del_positions_counts.at[n, 'count']
                    df_output_del_positions_counts.loc[n,'count'] = old + c1_list[i]

                if output_round == "r2":
                    old = df_output_del_positions_counts.at[n, 'count']
                    df_output_del_positions_counts.loc[n,'count'] = old + c2_list[i]
                
                #get positions for histogram input
                old = df_input_del_positions_counts.at[n, 'count']
                df_input_del_positions_counts.loc[n,'count'] = old + c_input_list[i]
                break
        
input_del_positions_counts = df_input_del_positions_counts['count'].tolist()
output_del_positions_counts = df_output_del_positions_counts['count'].tolist()
input_del_positions = df_input_del_positions_counts['position'].tolist()

#safe as .csv files
with open(f'deletions_position_count_input_{output_round}.csv', 'w', newline='') as file:
    writer = csv.writer(file, delimiter=",")
    writer.writerow(["position", "count"])

with open(f'deletions_position_count_output_{output_round}.csv', 'w', newline='') as file:
    writer = csv.writer(file, delimiter=",")
    writer.writerow(["position", "count"])

for i in range(0, 189):
    with open(f'deletions_position_count_input_{output_round}.csv', 'a', newline='') as file:
        writer = csv.writer(file, delimiter=",")
        writer.writerow([f"{input_del_positions[i]}", f"{input_del_positions_counts[i]}"])
    with open(f'deletions_position_count_output_{output_round}.csv', 'a', newline='') as file:
        writer = csv.writer(file, delimiter=",")
        writer.writerow([f"{input_del_positions[i]}", f"{output_del_positions_counts[i]}"])


#_____________________________________________________________________________________________________________________________
#plot histograms with results from analysis

input_counts_raw = []
with open (f"deletions_position_count_input_{output_round}.csv", 'rt') as f:
    reader = csv.reader(f, delimiter=',')
    for row in reader: 
        input_counts_raw.append(row)
input_counts_raw.pop(0)

output_counts_raw = []
with open (f"deletions_position_count_output_{output_round}.csv", 'rt') as f:
    reader = csv.reader(f, delimiter=',')
    for row in reader: 
        output_counts_raw.append(row)
output_counts_raw.pop(0)

#+54 to rescale from the start
positions = [(int(x[0]) + 54) for x in input_counts_raw]
counts_input = [int(x[1]) for x in input_counts_raw]
counts_output = [int(x[1]) for x in output_counts_raw]

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'
   
print(color.BLUE + "Position analysis" + color.END)
print("Plotting histograms...")

#calculate frequencies for histogram and define input variables (can be used to get frequencies for distinct positions) 
input_del_frequencies = [(i / Nt_input)*100 for i in counts_input]
output_del_frequencies = [(i / Nt)*100 for i in counts_output]
#print(output_del_frequencies)
ticks = np.arange(50,250,10)
stop_color = '#4575b4'
round_name = f"{output_round}"

plt.figure(figsize=(8.3, 2.4))
plt.bar(positions, input_del_frequencies, color='dimgray', linewidth=0, alpha=0.5, label="input library", width=1)
plt.bar(positions, output_del_frequencies, color=f'{stop_color}', linewidth=0, alpha=0.9, label=f"output {round_name}", width = 0.5)
plt.legend(loc="upper right")
plt.title(f"Single nucleotide deletions: positional frequency distribution after {round_name}", fontsize = 10)
plt.xlabel("Position", fontsize = 8)
plt.ylabel("Deletion frequnecy [%]", fontsize = 8)
plt.xticks(ticks, fontsize=8)
plt.yticks(fontsize=8)
plt.savefig(f"stop_position_bar_inputvsoutput_total_{round_name}.pdf")
plt.close()





