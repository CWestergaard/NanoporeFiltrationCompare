"""
Program: CoverageTimesort
Description: Simulate real-time data generation by sorting input fastq.gz file 
             into smaller files, where each file contains enough reads to reach
             a desired coverage, based on genome-size. 
             Reads are sorted based on their timestamp (Nanopore)
Version: 1.0
Author: Casper Westergaard

Example string:
python3 CoverageTimesort.py -i /srv/data/AS/CASW/data/q8/CPO20160077/CPO20160077.chop.q8.fastq.gz
 -o /srv/data/AS/CASW/data/q8/CPO20160077/ -gs 5m -cl 1,2,3,4,5,6,7,8,9,10,15,20,30,40,50,60
"""

#Import libraries
import datetime
import gzip
import sys
import os
import argparse
import shutil

###########################################################################
# GET INPUT
###########################################################################

# Input from commandline
parser = argparse.ArgumentParser(
        description='Sort input fastq.gz file into smaller files, where each file'
                    ' contains enough reads to reach a desired coverage based on'
                    ' genome-size. Reads are sorted based on their timestamp')
parser.add_argument('-i', type=str, dest='input_filename', 
                    help='Input fastq.gz file containing timestamped reads', required=True)
parser.add_argument('-o', type=str, dest='output_folder', 
                    help='Path to output folder', required=True)
parser.add_argument('-gs', type=str, dest='genome_size_string', 
                    help='Size of genome', required=True)
parser.add_argument('-cl', type=str, dest='coverage_list', 
                    help='List of desired coverages, a new file will be created for each element,'
                    ' coverages should be comma-delimitered and ordered from lowest to highest', required=True)
args = parser.parse_args()

# Check that input arguments are valid
try:
    infile = gzip.open(args.input_filename, 'r')
except IOError as error:
    sys.exit('Can\'t open file, reason:',str(error),'\n') 
args.output_folder += '/'
os.makedirs('{}timesorted_coverages/'.format(args.output_folder), exist_ok=True)

# Convert input genome size to int
if args.genome_size_string[-1].upper() == 'K':
    genome_size = int(args.genome_size_string[0:-1]) * 1000
elif args.genome_size_string[-1].upper() == 'M':
    genome_size = int(args.genome_size_string[0:-1]) * 1000000
elif args.genome_size_string[-1].upper() == 'G':
    genome_size = int(args.genome_size_string[0:-1]) * 1000000000
else:
    genome_size = int(args.genome_size_string)
    
# Check thal all coverages in list are ints
try:
    args.coverage_list = [int(x) for x in args.coverage_list.split(',')]
except ValueError as err:
    print('Coverages must be in integer-values:', str(err))
    sys.exit(1)
for i, cov in enumerate(args.coverage_list):
    if i == len(args.coverage_list)-1:
        break
    else:
        if cov > args.coverage_list[i + 1]:
            sys.exit('Coverages must be ordered from lowest to highest')

###########################################################################
# Load all reads into memory and sort them based on their timestamp
###########################################################################  
              
# Go through each read, save the timestamp, number of bases and the total read
reads = list()
flag = 0
total_bases = 0
line_count = 4
for line in infile:
    # Sequence
    if flag == 1:
        bases = len(line)-1
        total_bases += bases
        reads[-1][1] = bases
        flag = 0
    # Header
    if line_count == 4:
        reads.append([None,None,b''])
        if line.find(b'sample') != -1:	#High Accuracy Basecalling changes the header, and winter/summertime
                index = 5
                hour_start = line.find(b'T')
                hour_end = line.find(b':')
                new_hour = str(int(line[hour_start+1:hour_end])+1).encode('ascii')
                line = line[:hour_start+1]+new_hour+line[hour_end:]
        else:
                index = 4        
        time = line.split()[index].split(b'=')[1].decode('ascii')   # Get time from timestamp in header
        reads[-1][0] = datetime.datetime.strptime(time, '%Y-%m-%dT%H:%M:%SZ')
        flag = 1
        line_count = 0
    # Save all information, for output later
    reads[-1][2] += line
    line_count += 1
infile.close()

# Sort reads based on time        
reads.sort()

# Find required amount of reads to reach wanted coverage
basepairs = 0
coverage_reads = list()
coverage_counter = 0       
previous_count = 0
for i in range(len(reads)):     
    if basepairs < args.coverage_list[coverage_counter]*genome_size:
        basepairs += reads[i][1]
    else:
        coverage_reads.append([args.coverage_list[coverage_counter],reads[previous_count:i],i])
        basepairs += reads[i][1]
        previous_count = i
        coverage_counter += 1
        if basepairs > args.coverage_list[-1]*genome_size:
            break       
        
###########################################################################
# Create output files
###########################################################################

# Create filenames
isolate_filename = os.path.basename(args.input_filename)
isolate = '.'.join(isolate_filename.split('.')[0:-2]) # To not include .fastq.gz
output_filenames = list()
for i in range(len(coverage_reads)):
    if coverage_reads[i][0] < 10:
        output_filenames.append('{}timesorted_coverages/{}.cov_00{}.fastq.gz'.format(args.output_folder,isolate,coverage_reads[i][0]))        
    elif coverage_reads[i][0] < 100:
        output_filenames.append('{}timesorted_coverages/{}.cov_0{}.fastq.gz'.format(args.output_folder,isolate,coverage_reads[i][0]))  
    else:
        output_filenames.append('{}timesorted_coverages/{}.cov_{}.fastq.gz'.format(args.output_folder,isolate,coverage_reads[i][0]))  

# Write to first file
outfile = gzip.open(output_filenames[0], 'w')
for read in coverage_reads[0][1]:
    outfile.write(read[2])
outfile.close()
for i in range(1,len(coverage_reads)):
    # Copy last file and append to it
    shutil.copyfile(output_filenames[i-1], output_filenames[i])
    with gzip.open(output_filenames[i], 'a') as outfile:
        for read in coverage_reads[i][1]:
            outfile.write(read[2])


# Write coverage time information file
outfilename = ('{}timesorted_coverages/{}.cov_stats.txt'
                .format(args.output_folder, isolate))
with open(outfilename, 'w') as outfile:
    outfile.write('Input file: ' + args.input_filename)
    outfile.write('\nTotal bases in file: {}'.format(total_bases))
    outfile.write('\nTotal bases coverage: {}'.format(total_bases/genome_size))
    outfile.write('\nGenome size: {}'.format(args.genome_size_string))
    outfile.write('\nTime to reach coverage (Hours:Minutes:Seconds):\n')
    for i in range(len(coverage_reads)):
        if args.coverage_list[i] < 100:
            outfile.write('0')
        if args.coverage_list[i] < 10:
            outfile.write('0')
        time_since_start = (datetime.datetime.min + (reads[(coverage_reads[i][2])-1][0]-reads[0][0])).time()
        outfile.write('{}: {} \n'.format(args.coverage_list[i],time_since_start))
    