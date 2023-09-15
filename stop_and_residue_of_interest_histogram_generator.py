#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import SeqIO #can read and write most sequence formats
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

#_____________________________________________________________________________________________________________________________
parser = argparse.ArgumentParser(description=""" --- """)

parser.add_argument('-T', '--threads', type=int, default=0, help='Number of threads to execute in parallel. Defaults to CPU count.')
parser.add_argument('--variant_merge_file', help='output from dimsum with correct amino acid sequence and variant counts', required=True)
parser.add_argument('--round', help='enter selection round: either "r1" or "r2"', required=True)
parser.add_argument('--residue', help='enter residue of interest here - 1 letter code e.g "C"', required=True)

#extension for parser and threads
args = parser.parse_args()
threads = args.threads
if threads == 0:
    threads = multiprocessing.cpu_count()

#_____________________________________________________________________________________________________________________________
#set up input variables
variant_merge_file = args.variant_merge_file
output_round = str(args.round)
residue = str(args.residue)

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
sequences = [str(x[1]) for x in variant_merge_list]
c_input_list = [int(x[8]) for x in variant_merge_list]
c1_list = [int(x[9]) for x in variant_merge_list]
c2_list = [int(x[10]) for x in variant_merge_list]

#_____________________________________________________________________________________________________________________________
#Nt = total number of reads for a certain sequence
Nt_input = 0
Nt_1 = 0
Nt_2 = 0

for i in range(0, len(c_input_list)):
    Nt_input += c_input_list[i]
    Nt_1 += c1_list[i]
    Nt_2 += c2_list[i]
    #Nt_3 += c3_list[i]

if output_round == "r1":
    Nt = Nt_1
if output_round == "r2":
    Nt = Nt_2
#_____________________________________________________________________________________________________________________________
#generate data frames for histograms
data_frame_positions = list(np.arange(19, 84, 1))
data_frame_initial_fill = []
for i in range(0, 65):
    data_frame_initial_fill.append(0)

df_input_preliminary_stop_positions_counts = pd.DataFrame(list(zip(data_frame_positions, data_frame_initial_fill)), columns = ['position', 'count'])
df_input_residue_positions_counts = pd.DataFrame(list(zip(data_frame_positions, data_frame_initial_fill)), columns = ['position', 'count'])
df_output_preliminary_stop_positions_counts = pd.DataFrame(list(zip(data_frame_positions, data_frame_initial_fill)), columns = ['position', 'count'])
df_output_residue_positions_counts = pd.DataFrame(list(zip(data_frame_positions, data_frame_initial_fill)), columns = ['position', 'count'])

for i in range(0, len(sequences)):
    preliminary_stop = "F"
    residue_status = "F"
    sequence = sequences[i]
    for n in range(len(sequence)):
        #if there is a preliminary stop codon (before position 64, first 18 nucs are not sequenced)
        if sequence[n] == "*" and n < 65:
            if output_round == "r1":
                #get old value in dataframe, add current count and then replace the old with the new one
                old = df_output_preliminary_stop_positions_counts.at[n, 'count']
                df_output_preliminary_stop_positions_counts.loc[n,'count'] = old + c1_list[i]

            if output_round == "r2":
                old = df_output_preliminary_stop_positions_counts.at[n, 'count']
                df_output_preliminary_stop_positions_counts.loc[n,'count'] = old + c2_list[i]
                
            #get positions for histogram input
            old = df_input_preliminary_stop_positions_counts.at[n, 'count']
            df_input_preliminary_stop_positions_counts.loc[n,'count'] = old + c_input_list[i]
            
            #leave loop if preliminary stop was found -> double stops are igonored
            break
        
    #look if there are residues of interest 
    for n in range(len(sequence)):
        
        #if there is at least one residue of interest residue in the coding sequence of normal length
        if sequence[n] == f"{residue}" and n < 65:

            #generate data for input and output library histograms
            old = df_input_residue_positions_counts.at[n, 'count']
            df_input_residue_positions_counts.loc[n,'count'] = old + c_input_list[i]
            if output_round == "r1":
                old = df_output_residue_positions_counts.at[n, 'count']
                df_output_residue_positions_counts.loc[n,'count'] = old + c1_list[i]
            if output_round == "r2":
                old = df_output_residue_positions_counts.at[n, 'count']
                df_output_residue_positions_counts.loc[n,'count'] = old + c2_list[i]

        #only counts residues of interest in the coding sequence as it would break when there is a stop codon
        if sequence[n] == "*":
            break
#_____________________________________________________________________________________________________________________________
#print results from analysis

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

input_preliminary_stop_positions_counts = df_input_preliminary_stop_positions_counts['count'].tolist()
input_residue_positions_counts = df_input_residue_positions_counts['count'].tolist()
output_preliminary_stop_positions_counts = df_output_preliminary_stop_positions_counts['count'].tolist()
output_residue_positions_counts = df_output_residue_positions_counts['count'].tolist()

input_preliminary_stop_positions_frequencies = [(i / Nt_input)*100 for i in input_preliminary_stop_positions_counts]
input_residue_positions_frequencies = [(i / Nt_input)*100 for i in input_residue_positions_counts]
output_preliminary_stop_positions_frequencies = [(i/ Nt)*100 for i in output_preliminary_stop_positions_counts]
output_residue_positions_frequencies = [(i / Nt)*100 for i in output_residue_positions_counts]

positions = np.arange(19, 84, 1)
ticks = np.arange(15,90,5)


if output_round == "r1":
    residue_color = '#f46d43'
    stop_color = '#74add1'
    round_name = "round 1"
    
if output_round == "r2":
    residue_color = '#d73027'
    stop_color = '#4575b4'
    round_name = "round 2"
    
    
plt.figure(figsize=(9.4, 4))
plt.bar(positions, input_preliminary_stop_positions_frequencies, color='dimgray', linewidth=0, alpha=0.7, label="input library", width=1)
plt.bar(positions, output_preliminary_stop_positions_frequencies, color=f'{stop_color}', linewidth=0, alpha=0.9, label=f"output {round_name}", width = 0.5)
plt.legend(loc="upper right")
plt.title(f"Trunated variants: positional frequency distribution after {round_name}")
plt.xlabel("Position")
plt.ylabel("Fraction of total count number [%]")
plt.xticks(ticks)
plt.savefig(f"stop_position_bar_inputvsoutput_total_{round_name}.pdf")
plt.close()

plt.figure(figsize=(9.4, 4))
plt.bar(positions, input_residue_positions_frequencies, color='dimgray', linewidth=0, alpha=0.7, label="input library", width = 1)
plt.bar(positions, output_residue_positions_frequencies, color=f'{residue_color}', linewidth=0, alpha=0.9, label=f"output {round_name}", width = 0.5)
plt.legend(loc="upper right")
plt.title(f"{residue} residues: positional frequency distribution after {round_name}")
plt.xlabel("Position")
plt.ylabel("Fraction of total count number [%]")
plt.xticks(ticks)
plt.savefig(f"{residue}_position_normal_bar_inputvsoutput_total_{round_name}.pdf")
plt.tick_params(axis='x', which='minor')
plt.close()




