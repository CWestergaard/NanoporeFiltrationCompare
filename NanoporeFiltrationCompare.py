"""
Program: NanoporeFiltrationCompare
Description: Wrapper to run the full Nanopore filtration and AMR comparison with reference pipeline.
Version: 1.0
Author: Casper Westergaard
Example string:
python3 /srv/data/AS/CASW/scripts/Wrapper/GitHub/NanoporeFiltrationCompare.py -i /srv/data/AS/CASW/data/q8_subset/ -o /srv/data/AS/CASW/resultater/q8_subset/ -gs 5m -cl 1,2,3,4,5,6,7,8,9,10,15,20,30,40,50,60 -md 0.2 -df 0.45 -ip /srv/data/AS/CASW/data/phenotypes.txt -op /srv/data/AS/CASW/resultater/Plots -refdb /srv/data/DB/kma/resfinder_db -sct 1
"""

#Import Libraries
import os
import argparse
import glob
import subprocess

###########################################################################
# GET INPUT
###########################################################################
# Input from commandline
# Required input
parser = argparse.ArgumentParser(description='Wrapper to run the full Nanopore filtration and AMR comparison with reference pipeline.')
parser.add_argument('-ip', type=str, dest='input_phenotypes', 
                    help='Input file containing name and phenotype of all genes in database', required=True)
parser.add_argument('-o', type=str, dest='output_folder', 
                    help='Path to folder to put result in', required=True)
parser.add_argument('-i', type=str, dest='input_folder', 
                    help='Path to folder containing individual folders for all isolates', required=True)
parser.add_argument('-op', type=str, dest="output_plots", 
                    help='Folder to output plots in', required=True)
parser.add_argument('-gs', type=str, dest="genome_size_string",
                    help='Size of genome', required=True)
parser.add_argument('-cl', type=str, dest="coverage_list", 
                    help='List of desired coverages, a new file will be created for each element,'
                    'coverages should be comma-delimitered and ordered from lowest to highest', required=True)
parser.add_argument('-refdb', type=str, dest='refdb',
                    help='Path to ResFinder database', required=True)
# Optional input
parser.add_argument('-df', type=float, dest='depth_filtration', default='0',
                    help='Set depth filtration, default 0, input value between 0 and 1', required=False)
parser.add_argument('-md', type=float, dest='min_depth', default='0', 
                    help='Minimum depth for gene-groups to be included', required=False)
parser.add_argument('-id', type=float, dest='identity_cutoff', default='90',
                    help='Minimum identity for genes to be included, default 90, input value between 0 and 100', required=False)
parser.add_argument('-cov', type=float, dest='coverage_cutoff', default='90',
                    help='Minimum coverage for genes to be included, default 90, input value between 0 and 100', required=False)
parser.add_argument('-sct', type=int, dest='skip_coverage_timesort', default='0', 
                    help='Skip coverage_timesort script, if files have already been created. -sct 1 to skip', required=False)
parser.add_argument('-sk', type=int, dest='skip_kma', default='0', 
                    help='Skip kma script, if kma has already been run on files. -sk 1 to skip', required=False)
args = parser.parse_args()

args.output_folder += '/'
args.input_folder += '/'
args.output_plots += '/'

os.chdir(args.input_folder)    
isolate_folders = glob.glob('*')
dirname = os.path.dirname(os.path.abspath(__file__))

###########################################################################
# CALL CoverageTimesort.py
###########################################################################

if args.skip_coverage_timesort == 0:
# Call coverage_timesort script for all isolates    
    print('STEP 1 - Coverage Timesort')
    for i in range(len(isolate_folders)):
        output_path = '{}{}'.format(args.input_folder,isolate_folders[i])
        input_file = '{}{}/{}.chop.q8.fastq.gz'.format(args.input_folder,isolate_folders[i],isolate_folders[i])
        if not os.path.isfile(input_file):
            input_file = '{}{}/{}.q8.fastq.gz'.format(args.input_folder,isolate_folders[i],isolate_folders[i]) # No .chop

        script_path = os.path.join(dirname, 'CoverageTimesort.py')   
        print('Running script on: {}'.format(isolate_folders[i]))
        cmd_CT = ['python3', script_path, '-i', input_file, '-o', output_path,
                  '-gs', args.genome_size_string, '-cl', args.coverage_list]
        subprocess.run(cmd_CT, check=True)
        print('Script finished with: {}'.format(isolate_folders[i]))
        
elif args.skip_coverage_timesort == 1:
    print('Skipping STEP 1 - Coverage Timesort')

###########################################################################
# Call RunKMA.py
###########################################################################

if args.skip_kma == 0:
# Call run_kma_on_folder script for all isolates
    print('STEP 2 - KMA')
    for i in range(len(isolate_folders)):
        output_path = '{}{}/'.format(args.output_folder,isolate_folders[i])  
        input_folder = '{}{}/'.format(args.input_folder,isolate_folders[i])

        script_path = os.path.join(dirname, 'RunKMA.py')
        print('Running script on: {}'.format(isolate_folders[i]))
        cmd_KMA = ['python3', script_path, '-i', input_folder, '-o', output_path ,'-refdb', args.refdb]
        subprocess.run(cmd_KMA, check=True)
        print('Script finished with: {}'.format(isolate_folders[i]))
        
elif args.skip_kma == 1:
    print('Skipping STEP 2 - KMA')

###########################################################################
# CALL NanoporeFiltration.py
###########################################################################

# Call Nanopore filtration script for all isolates
print('STEP 3 - Filtration')
for i in range(len(isolate_folders)):
    output_path = '{}{}/'.format(args.output_folder,isolate_folders[i]) 
    input_path = '{}{}/kma/'.format(args.output_folder,isolate_folders[i])
    ref_file = glob.glob('{}{}/*.tsv'.format(args.input_folder, isolate_folders[i]))

    script_path = os.path.join(dirname, 'NanoporeFiltration.py')
    print('Running script on: {}'.format(isolate_folders[i]))
    cmd_NF = ['python3', script_path, '-ip', args.input_phenotypes, 
                '-if', input_path, '-o', output_path, '-ref', ref_file[0], '-df', 
                str(args.depth_filtration), '-md', str(args.min_depth), '-cl', args.coverage_list]
    subprocess.run(cmd_NF, check=True)
    print('Script finished with: {}'.format(isolate_folders[i]))

###########################################################################
# CALL ComparisonPlots.py
###########################################################################
    
# Call plot script and craete plots for all isolates
print('STEP 4 - Plots')
script_path = os.path.join(dirname, 'ComparisonPlots.py')
cmd_plot = ['python3', script_path, '-i', args.output_folder, '-o', args.output_plots,
            '-df', str(args.depth_filtration), '-md', str(args.min_depth)]
subprocess.run(cmd_plot, check=True)
print('Pipeline done')


    
        

