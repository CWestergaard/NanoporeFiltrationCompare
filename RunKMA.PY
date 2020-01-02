#!/usr/bin/env python3
"""
Program: RunKMA
Description: Runs kma on all fastq.gz files in the given input folder, using the ResFinder database.
             Illumina paired end reads should be in the input folder. Nanopore coverage files
             should be in a subfolder named timesorted_coverages.
Version: 1.0
Author: Casper Westergaard

Example string:
python3 /srv/data/AS/CASW/scripts/Wrapper/GitHub/RunKMA.py -i /srv/data/AS/CASW/data/HAC/CPO20180005 -o /srv/data/AS/CASW/resultater/HAC/CPO20180005 -refdb /srv/data/DB/kma/resfinder_db
"""

#Import libraries
import os
import sys
import subprocess
import glob
import argparse

###########################################################################
# GET INPUT
###########################################################################

# Input from commandline
parser = argparse.ArgumentParser(description='Runs kma with the ResFinder database on all fastq.gz files in the given input folder')
parser.add_argument('-i', type=str, dest='input_folder', 
                    help='Isolate folder containing the fastq.gz files to run kma on', required=True)
parser.add_argument('-o', type=str, dest='output_folder', 
                    help='Output folder to put kma results in', required=True)
parser.add_argument('-refdb', type=str, dest='refdb',
                    help='Path to ResFinder database', required=True)
args = parser.parse_args()

db = 'resfinder'
#refdb = '/srv/data/DB/kma/resfinder_db'

# Check if input folder exists
if not os.path.exists(args.input_folder):
    print('Input folder does not exist')
    sys.exit(1)

# Check if timesorted_coverages exists in input folder
if not os.path.exists(args.input_folder+'/timesorted_coverages'):
    print('The folder timesorted_coverages can not be found in the input folder')
    sys.exit(1)

# Create output folder
os.makedirs('{}{}'.format(args.output_folder,'kma/'), exist_ok=True)

###########################################################################
# RUN KMA
###########################################################################  
      
# Run on Illumina paired end files
os.chdir(args.input_folder)
r1 = ''
r2 = ''
for file in glob.glob('*.fastq.gz'):
    if 'R1' in file:
        r1 = file
    elif 'R2' in file:
        r2 = file
if r1 != '' and r2 != '':
    cmd_kma = ['kma', '-ipe', args.input_folder+r1,
               args.input_folder+r2, '-o',
               args.output_folder+'kma/kma_'+db+'.'+r1.split('.')[0].split('_')[0]+'.z.illumina.trimmed',
               '-t_db', args.refdb, '-1t1']
    subprocess.run(cmd_kma, check=True)

# Run on Nanopore timesorted coverage files
os.chdir(args.input_folder+'timesorted_coverages')
for file in glob.glob('*.fastq.gz'):
    full_path_file = args.input_folder+'timesorted_coverages/'+file
    cmd_kma = ['kma', '-i', full_path_file, '-o', args.output_folder+'kma/kma_'+db+'.'+file,
               '-t_db', args.refdb, '-mem_mode', '-mp', '20', '-bcNano', '-t', '4']
    subprocess.run(cmd_kma, check=True)
