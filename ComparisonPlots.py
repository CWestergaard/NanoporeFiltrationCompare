"""
Program: ComparisonPlots
Description: Script for generating comparison plots for nanopore files of varied coverages, compared to a reference.
             Input is based on results from the script NanoporeFiltration.py
             Creates plots for antibiotic class, phenotype, genes and gene-group for each inputfile, as well as summary plots.
Version: 1.0
Author: Casper Westergaard
Example string:
python3 /srv/data/AS/CASW/scripts/Wrapper/GitHub/ComparisonPlots.py -i /srv/data/AS/CASW/resultater/q8_subset/ -o /srv/data/AS/CASW/resultater/Plots/ -df 0.45 -md 0.2
"""

#Import libraries
import argparse
import glob
import sys
import os
import matplotlib.pyplot as plt
plt.switch_backend('agg')

###########################################################################
# FUNCTIONS
########################################################################### 

def count_results(labels, results,correct,missing,extra,extra_and_missing):
    """Count correct, missing, extra and extra_and_missing for each label in given input dictionary"""
    for label in labels:
        correct.append(0),missing.append(0),extra.append(0),extra_and_missing.append(0)
        for result in results[label]:
            if result[0] > 0 and result[1] == 0 and result[2] == 0:
                correct[-1] += 1
            elif result[1] > 0 and result[2] == 0:
                missing[-1] += 1
            elif result[2] > 0 and result[1] == 0:
                extra[-1] += 1
            else:
                extra_and_missing[-1] += 1
    return correct,missing,extra,extra_and_missing

def plot(isolate_names, labels, results, title, folder_name):
    """ Create plots of Correct, Extra and Missing hits for each isolate.
        Also create summary plot for all isolates"""
    #Create individual plots for each isolate
    pos = list(range(len(labels)))
    width = 0.25
    for i in range(len(isolate_names)):
        fig, ax = plt.subplots(figsize=(10,5))
        correct = [results[label][i][0] for label in labels]
        missing = [results[label][i][1] for label in labels]
        extra = [results[label][i][2] for label in labels]
        bar1 = plt.bar(pos, correct, width, alpha=0.5, color='#28a122')
        bar2 = plt.bar([p + width for p in pos], missing, width, alpha=0.5, color='#3300FF')
        bar3 = plt.bar([p + width*2 for p in pos], extra, width, alpha=0.5, color='#FFC222')
        ax.set_ylabel('Number of ' + title)
        ax.set_xlabel('Testfiles')
        ax.set_title('{} - Comparison of {}  between reference and testfiles'.format(isolate_names[i],title))
        ax.set_xticks([p + 1.5 * width for p in pos])
        plt.xticks(rotation=270)
        ax.set_xticklabels(labels)
        plt.xlim(min(pos)-width, max(pos)+width*4)
        plt.ylim([0, max(correct + missing + extra)+10])
        for rect in bar1 + bar2 + bar3:
            height = rect.get_height()
            plt.text(rect.get_x() + rect.get_width()/2.0, height, '%d' % int(height), ha='center', va='bottom')
        plt.legend(['Correct hits', 'Missing hits', 'Extra hits'], loc='upper left')
        plt.figtext(0.92, 0.5, filtrations[i])
        plt.grid()
        plt.savefig('{}/{}_{}.png'.format(output_path,isolate_names[i],title),bbox_inches='tight')
        plt.close(fig)
        
    #Create summary plot of all isolates
    correct, missing, extra, extra_and_missing = [list(),list(),list(),list()]
    count_results(labels, results,correct,missing,extra,extra_and_missing)
    pos = [i * 1.5 for i in pos]
    fig, ax = plt.subplots(figsize=(10,5))
    bar1 = plt.bar(pos, correct, width, alpha=0.5, color='#28a122')
    bar2 = plt.bar([p + width for p in pos], missing, width, alpha=0.5, color='#3300FF')
    bar3 = plt.bar([p + width*2 for p in pos], extra, width, alpha=0.5, color='#FFC222')               
    bar4 = plt.bar([p + width*3 for p in pos], extra_and_missing, width, alpha=0.5, color='#994C00')
    ax.set_ylabel('Number of Isolates')
    ax.set_xlabel('Testfiles')
    ax.set_title('Summary of results from all isolates - {}'.format(title))
    ax.set_xticks([p + 1.5 * width for p in pos])
    plt.xticks(rotation=270)
    ax.set_xticklabels(labels)
    plt.xlim(min(pos)-width, max(pos)+width*4)
    plt.ylim([0, max(correct + missing + extra + extra_and_missing)+10])
    for rect in bar1 + bar2 + bar3 + bar4:
        height = rect.get_height()
        if height > 0:
            plt.text(rect.get_x() + rect.get_width()/2.0, height, '%d' % int(height), ha='center', va='bottom')
    plt.legend(['Only Correct hits', 'Correct + Missing hits', 'Correct + Extra hits', 'Correct + Missing + Extra Hits'], loc='upper left')
    plt.figtext(0.92, 0.5, filtrations[0])
    plt.grid()
    file_name = '{}_{}_Multiple_Isolates.png'.format(folder_name,title)
    plt.savefig('{}/{}'.format(output_path,file_name),bbox_inches='tight')
    plt.close(fig)
    
###########################################################################
# GET INPUT
###########################################################################
    
# Input from commandline
# Required input
parser = argparse.ArgumentParser(description='Create plots based on information in summary_stats file for each isolate')
parser.add_argument('-i', type=str, dest='input_folder',
                    help='Path to input folder, each isolate should have its own folder', required=True)
parser.add_argument('-o', type=str, dest='output',
                    help="Path to folder to put plots in", required=True)
parser.add_argument('-df', type=float, dest='depth_filtration', default='0',
                    help='Depth filtration used for generating input', required=False)
parser.add_argument('-md', type=float, dest='min_depth', default='0',
                    help='Minimum depth for gene-groups to be included', required=False)
args = parser.parse_args()

input_string = '{}**/Depth={}_Min_Depth={}/*_summary_stats.txt'.format(args.input_folder,args.depth_filtration,args.min_depth)

# Depending if input_folder has one or two /
folder_name = args.input_folder.split('/')[-2]
if folder_name == '':
    folder_name = args.input_folder.split('/')[-3]

output_path = '{}Depth={}_Min_Depth={}/'.format(args.output,args.depth_filtration,args.min_depth)    
os.makedirs('{}'.format(output_path), exist_ok=True)
 
###########################################################################
# LOAD INPUT DATA
###########################################################################
     
# Load results from isolates summary file
filtrations = list()
results_class = dict()
results_phenotype = dict()
results_gene_group = dict()
results_gene = dict()  
isolate_names = list()
for file in glob.glob(input_string, recursive=True):
    isolate_names.append(file.split('/')[-1].split('_')[0]) # Get isolate name based on filename
    filtrations.append('')
    #input_file = open(file,'r')
    with open(file, 'r') as input_file:
        filtration_flag, class_flag, phenotype_flag, gene_group_flag, gene_flag, labels = [1,0,0,0,0,list()]
        for line in input_file:
            if line == '\n':
                filtration_flag, class_flag, phenotype_flag, gene_group_flag, gene_flag, labels = [0,0,0,0,0,list()] 
            elif line.startswith('##'):
                                 header = line
            elif line == '#Antibiotic class\n':
                class_flag = 1
            elif line == '#Phenotype\n':
                phenotype_flag = 1   
            elif line == '#Gene-group\n':
                gene_group_flag = 1  
            elif line == '#Gene\n':
                gene_flag = 1              
            elif filtration_flag == 1:
                filtrations[-1] += line
            else:
                line = line.split('\t')
                label = line[0]
                correct_hit = int(line[1])
                missing_hit = int(line[2])
                extra_hit = int(line[3])
                labels.append(label)
                if class_flag == 1:
                    results_class.setdefault(label, []).append([correct_hit,missing_hit,extra_hit])               
                elif phenotype_flag == 1:
                    results_phenotype.setdefault(label, []).append([correct_hit,missing_hit,extra_hit])
                elif gene_group_flag == 1:
                    results_gene_group.setdefault(label, []).append([correct_hit,missing_hit,extra_hit])
                elif gene_flag == 1:
                    results_gene.setdefault(label, []).append([correct_hit,missing_hit,extra_hit])

#Check that filtrations and labels are identical for all isolates                    
if len(set(filtrations)) > 1:
    print('The filtrations are not equal for all input files.')
    sys.exit(1)
for label in results_class.keys():
    if len(results_class[label]) != len(glob.glob(input_string, recursive=True)):
        print('Labels are not the same for all input files.')
        print(label)
        sys.exit(1)

###########################################################################
# CALL FUNCTIONS TO GENERATE PLOTS
###########################################################################

#Call function plot, to create plots for individual isolates aswell as summary plots
plot(isolate_names,labels,results_class,'Antimicrobial classes', folder_name)
plot(isolate_names,labels,results_phenotype,'Phenotypes', folder_name)
plot(isolate_names,labels,results_gene_group,'Gene-groups', folder_name)
plot(isolate_names,labels,results_gene,'Gene-variants', folder_name)
