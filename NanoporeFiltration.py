"""
Program: NanoporeFiltration
Description: Apply filtration and compare kma results from testfiles (Reduced coverage) with reference results from ABRicate.
             Output effect of filtration, results after filtration and count of correct, missing and extra hits for plots.
Version: 1.0
Author: Casper Westergaard
Example string:
python3 /srv/data/AS/CASW/scripts/Wrapper/GitHub/NanoporeFiltration.py -ip /srv/data/AS/CASW/data/phenotypes.txt -if /srv/data/AS/CASW/resultater/HAC/CPO20180005/kma -o /srv/data/AS/CASW/resultater/HAC/CPO20180005/ -ref /srv/data/AS/CASW/data/HAC/CPO20180005/assembly.fasta.resfinder.tsv -df 0.45 -md 0.2 -cl 1,2,3,4,5,6,7,8,9,10,15,20,30,40,50,60
"""
#Import libraries
import argparse
import csv
import glob
import os
import sys

###########################################################################
# FUNCTIONS
########################################################################### 

class Gene:
    def __init__(self, ab_class, phenotype, group):
        self.ab_class = ab_class
        self.phenotype = phenotype
        self.group = group
        
def get_genes_from_ref_file(csvfile, coverage_cutoff, identity_cutoff):
    """ Object stores a list of genes that is above the given identity and
        coverage thresholds.
        Each gene is a string created from gene_name and accession no.
    """
    genes = list()

    with open(csvfile, 'r') as ref_results_file:

        reader = csv.reader(ref_results_file, delimiter='\t')
        # ignore first row
        next(reader)

        for row in reader:
            coverage = float(row[8])
            identity = float(row[9])
            gene_name = row[4]
            accesion = row[11]

            if(coverage > coverage_cutoff and identity > identity_cutoff):
                genes.append("{}_{}".format(gene_name, accesion))

    return genes

def get_genes_from_test_file(csvfile, coverage_cutoff, identity_cutoff):
    """ Object stores a list of genes that is above the given identity and
        coverage thresholds.
        Each gene is a string created from gene_name and accession no.
    """
    genes = list()

    with open(csvfile, 'r') as test_results_file:

        reader = csv.reader(test_results_file, delimiter='\t')
        # ignore first row
        next(reader)

        for row in reader:
            template_identity = float(row[4])
            template_coverage = float(row[5])
            query_identity = float(row[6])
            query_coverage = float(row[7])
            depth = float(row[8])
            gene_name = row[0].split('~~~')[1]
            accesion = row[0].split('~~~')[2]

            if (template_identity >= identity_cutoff and template_coverage >= coverage_cutoff \
            and query_identity >= identity_cutoff and query_coverage >= coverage_cutoff):
                genes.append(["{}_{}".format(gene_name, accesion),[template_identity,template_coverage,query_identity,query_coverage,depth]])

    return genes
        
def get_gene_group(gene_name):
    """ Extract the gene group from the gene name and returns it
    """
    gene_group = gene_name[0:3]

    if gene_group == 'bla':
        if len(gene_name.split('-')) == 1:

            if gene_name[0:6] == 'blaBEL':
                gene_group = 'blaBEL'
            else:
                gene_group = gene_name.split('_')[0]

        else:
            gene_group = gene_name.split('-')[0]

    return gene_group


def group_genes(gene_name, gene_group_dict):
    """ Takes a gene name and adds it to the appropriate gene group in the
        given dictionary
    """
    gene_group = get_gene_group(gene_name)

    gene_list = gene_group_dict.get(gene_group, [])
    gene_list.append(gene_name)
    gene_group_dict[gene_group] = gene_list

    return gene_group_dict

###########################################################################
# GET INPUT
###########################################################################
    
# Input from commandline
# Required input
parser = argparse.ArgumentParser(description='Compare results from reference and testfiles')
parser.add_argument('-ip', type=str, dest='input_phenotypes',
                    help='Input file containing name and phenotype of all genes in database', required=True)
parser.add_argument('-o', type=str, dest='output',
                    help='Path to folder to put result in', required=True)
parser.add_argument('-ref', type=str, dest='ref_file',
                    help='Path to results from reference, used to compare other files with', required=True)
parser.add_argument('-if', type=str, dest='input_folder',
                    help='Path to folder containing kma results files for comparison with reference results', required=True)
parser.add_argument('-cl', type=str, dest='coverage_list',
                    help='List of desired coverages, a new file will be created for each element,'
                    ' coverages should be comma-delimitered and ordered from lowest to highest', required=True)

#Optional input
parser.add_argument('-df', type=float, dest='depth_filtration', default='0',
                    help='Set depth filtration, default 0, input value between 0 and 1', required=False)
parser.add_argument('-md', type=float, dest='min_depth', default='0',
                    help='Minimum depth for gene-groups to be included, default 0, input value between 0 and 1', required=False)
parser.add_argument('-id', type=float, dest='identity_cutoff', default='90',
                    help='Minimum identity for genes to be included, default 90, input value between 0 and 100', required=False)
parser.add_argument('-cov', type=float, dest='coverage_cutoff', default='90',
                    help='Minimum coverage for genes to be included, default 90, input value between 0 and 100', required=False)
args = parser.parse_args()

# Initialize variables
output_path = '{}Depth={}_Min_Depth={}/'.format(args.output,args.depth_filtration,args.min_depth)       
os.makedirs('{}'.format(output_path), exist_ok=True)
isolate = args.input_folder.split('/')[-3] # Get isolate name based on input folder

# Check thal all coverages in list are ints
try:
    args.coverage_list = [int(x) for x in args.coverage_list.split(',')]
except ValueError as err:
    print('Coverages must be in integer-values:', str(err))
    sys.exit(1)
for i, cov in enumerate(args.coverage_list):
    if i < len(args.coverage_list)-1:
        if cov > args.coverage_list[i + 1]:
            sys.exit('Coverages must be ordered from lowest to highest')

# Create labels and mininum depth relative to coverage
labels = list()
min_depths = list()
for i in range(len(args.coverage_list)):
    labels.append('Nanopore Cov_'+str(args.coverage_list[i]))
    if args.min_depth > 0:
        min_depths.append(args.min_depth*int(args.coverage_list[i])) 
    else:
        min_depths.append(0) 
# Add Illumina
labels.append('Illumina')  

###########################################################################
# GET GENES FROM REFERENCE AND TESTFILES
###########################################################################

# Create gene-database containing Antibiotic class, Phenotype and Gene-group for each gene
# Based on phenotypes.txt from the ResFinder database
inputfile = open(args.input_phenotypes,'r')
gene_dict_db = dict()
inputfile.readline()    # To not include the header in the list
for line in inputfile:
    line = line.split('\t')
    gene_name = line[0].strip()
    phenotype = line[2].replace(' ','')
    antibiotic_class = line[1].replace(' ','').lower()
    # Inconsistencies between phenotypes.txt and antibiotic_classes.txt 
    if antibiotic_class == 'lincosamides':
        antibiotic_class = 'lincosamide' 
    gene_group = get_gene_group(gene_name)
    gene_dict_db[gene_name] = Gene(antibiotic_class, phenotype, gene_group)

# Create database of genes in reference file
ref_genes = get_genes_from_ref_file(args.ref_file, args.coverage_cutoff, args.identity_cutoff)
ref_db = dict()
for ref_gene_name in ref_genes:
    if(ref_gene_name in gene_dict_db):
        gene = gene_dict_db[ref_gene_name]
        if gene.ab_class not in ref_db:
            ref_db[gene.ab_class] = dict()
        if gene.phenotype not in ref_db[gene.ab_class]:
            ref_db[gene.ab_class][gene.phenotype] = dict()
        if gene.group not in ref_db[gene.ab_class][gene.phenotype]:
            ref_db[gene.ab_class][gene.phenotype][gene.group] = list()
        ref_db[gene.ab_class][gene.phenotype][gene.group].append(ref_gene_name)
        
# Create database of genes in test files. Save depth information about each gene-group.
# Gene-group information is used for filtration later on.
testfiles_db = list()   # list of dictionaries, one for each file 
testfiles_gene_groups = list() # list of dictionaries, one for each file 
input_files = glob.glob(args.input_folder+'*.res')
input_files.sort()
for file in input_files:
    testfiles_db.append(dict())
    testfiles_gene_groups.append(dict())
    testfile_genes = get_genes_from_test_file(file, args.coverage_cutoff, args.identity_cutoff)
    for testfile_gene_name in testfile_genes:
        if(testfile_gene_name[0] in gene_dict_db):
            gene = gene_dict_db[testfile_gene_name[0]]
            if gene.ab_class not in testfiles_db[-1]:
                testfiles_db[-1][gene.ab_class] = dict()
            if gene.phenotype not in testfiles_db[-1][gene.ab_class]:
               testfiles_db[-1][gene.ab_class][gene.phenotype] = dict()
            if gene.group not in testfiles_db[-1][gene.ab_class][gene.phenotype]:
                testfiles_db[-1][gene.ab_class][gene.phenotype][gene.group] = list()
            testfiles_db[-1][gene.ab_class][gene.phenotype][gene.group].append(testfile_gene_name)
            if gene.group not in testfiles_gene_groups[-1]:
                testfiles_gene_groups[-1][gene.group] = list()
            testfiles_gene_groups[-1][gene.group].append(testfile_gene_name)

###########################################################################
# DEPTH FILTRATION ON NANOPORE DATA
########################################################################### 
            
# Sort genes in testfiles_gene_groups based on depth, and save total depth of each gene-group in total_depths
total_depths = list() 
for testfile in testfiles_gene_groups:
    total_depths.append(dict())
    for group in testfile:
        testfile[group].sort(key=lambda x: x[1][4],reverse = True)    # Sort genes based on depth
        total_depths[-1][group] = sum([testfile[group][x][1][4] for x in range(len(testfile[group]))])

# Filter by depth. If a gene has a depth lower than depth_cutoff*highest depth gene of that gene group,
# or a gene-group depth lower than min_depth, then discard it.
testfiles_db_filtered = list()
for i in range(len(testfiles_db)-1): # No filtration on Illumina
    testfiles_db_filtered.append(dict())
    for ab_class in testfiles_db[i]:
        for phenotype in testfiles_db[i][ab_class]:
            for group in testfiles_db[i][ab_class][phenotype]:
                for gene in testfiles_db[i][ab_class][phenotype][group]:
                    gene_depth = gene[1][4]
                    max_gene_depth = testfiles_gene_groups[i][group][0][1][4]  # Depth of highest depth gene in that gene-group
                    # Depth filtration
                    if gene_depth >= args.depth_filtration*max_gene_depth and total_depths[i][group] > min_depths[i]:
                        if ab_class not in testfiles_db_filtered[i]:
                            testfiles_db_filtered[i][ab_class] = dict()
                        if phenotype not in testfiles_db_filtered[i][ab_class]:
                            testfiles_db_filtered[i][ab_class][phenotype] = dict()
                        if group not in testfiles_db_filtered[i][ab_class][phenotype]:
                            testfiles_db_filtered[i][ab_class][phenotype][group] = list()
                        testfiles_db_filtered[i][ab_class][phenotype][group].append(gene)
testfiles_db_filtered.append(testfiles_db[-1])  # Add Illumina results to filtered results

###########################################################################
# COMPARE FILTRATED RESULTS WITH REFERENCE
########################################################################### 
 
# Compare testfiles with reference and save Correct, Missing and Extra hits.
# Planned update: Include functions to improve readability
missing_class = list()
extra_class = list()
correct_class = list()
missing_phenotype = list()
extra_phenotype = list()
correct_phenotype = list()
missing_group = list()
extra_group = list()
correct_group = list()
missing_gene = list()
extra_gene = list()
correct_gene = list()
for i in range(len(testfiles_db_filtered)):
    missing_class.append(dict())
    extra_class.append(dict())
    correct_class.append(dict())
    missing_phenotype.append(dict())
    extra_phenotype.append(dict())
    correct_phenotype.append(dict())
    missing_group.append(dict())
    extra_group.append(dict())
    correct_group.append(dict())
    missing_gene.append(dict())
    extra_gene.append(dict())
    correct_gene.append(dict())
    # Compare testfiles with reference to identify correct and extra hits
    for ab_class in testfiles_db_filtered[i]:
        for phenotype in testfiles_db_filtered[i][ab_class]:
            for group in testfiles_db_filtered[i][ab_class][phenotype]:
                for gene_info in testfiles_db_filtered[i][ab_class][phenotype][group]:  # Also contains information about depth, coverage and identity
                    gene = gene_info[0]     # Gene name and accesion number, to compare with reference.
                  
                    if ab_class in ref_db:
                        correct_class[i].setdefault(ab_class, []).append(gene_info)
                        
                        if phenotype in ref_db[ab_class]:
                            correct_phenotype[i].setdefault(ab_class, {})
                            correct_phenotype[i][ab_class].setdefault(phenotype, []).append(gene_info)
                            
                            if group in ref_db[ab_class][phenotype]:
                                correct_group[i].setdefault(ab_class, {})
                                correct_group[i][ab_class].setdefault(phenotype, {})
                                correct_group[i][ab_class][phenotype].setdefault(group, []).append(gene_info)
                                
                                if gene in ref_db[ab_class][phenotype][group]:
                                    correct_gene[i].setdefault(ab_class, {})
                                    correct_gene[i][ab_class].setdefault(phenotype, {})
                                    correct_gene[i][ab_class][phenotype].setdefault(group, {})
                                    correct_gene[i][ab_class][phenotype][group].setdefault(gene, []).append(gene_info)
                                else:
                                    extra_gene[i].setdefault(ab_class, {})
                                    extra_gene[i][ab_class].setdefault(phenotype, {})
                                    extra_gene[i][ab_class][phenotype].setdefault(group, {})
                                    extra_gene[i][ab_class][phenotype][group].setdefault(gene, []).append(gene_info)
                                    
                            else:
                                extra_group[i].setdefault(ab_class, {})
                                extra_group[i][ab_class].setdefault(phenotype, {})
                                extra_group[i][ab_class][phenotype].setdefault(group, []).append(gene_info)
                                extra_gene[i].setdefault(ab_class, {})
                                extra_gene[i][ab_class].setdefault(phenotype, {})
                                extra_gene[i][ab_class][phenotype].setdefault(group, {})
                                extra_gene[i][ab_class][phenotype][group].setdefault(gene, []).append(gene_info)
                                
                        else:
                            extra_phenotype[i].setdefault(ab_class, {})
                            extra_phenotype[i][ab_class].setdefault(phenotype, []).append(gene_info)
                            extra_group[i].setdefault(ab_class, {})
                            extra_group[i][ab_class].setdefault(phenotype, {})
                            extra_group[i][ab_class][phenotype].setdefault(group, []).append(gene_info)
                            extra_gene[i].setdefault(ab_class, {})
                            extra_gene[i][ab_class].setdefault(phenotype, {})
                            extra_gene[i][ab_class][phenotype].setdefault(group, {})
                            extra_gene[i][ab_class][phenotype][group].setdefault(gene, []).append(gene_info)
                            
                    else:
                        extra_class[i].setdefault(ab_class, []).append(gene_info)
                        extra_phenotype[i].setdefault(ab_class, {})
                        extra_phenotype[i][ab_class].setdefault(phenotype, []).append(gene_info)
                        extra_group[i].setdefault(ab_class, {})
                        extra_group[i][ab_class].setdefault(phenotype, {})
                        extra_group[i][ab_class][phenotype].setdefault(group, []).append(gene_info)
                        extra_gene[i].setdefault(ab_class, {})
                        extra_gene[i][ab_class].setdefault(phenotype, {})
                        extra_gene[i][ab_class][phenotype].setdefault(group, {})
                        extra_gene[i][ab_class][phenotype][group].setdefault(gene, []).append(gene_info)
                        
    # Compare reference with testfiles to identify missing hits                    
    for ab_class in ref_db:                        
        for phenotype in ref_db[ab_class]:
            for group in ref_db[ab_class][phenotype]:
                for gene in ref_db[ab_class][phenotype][group]:
                    
                    if ab_class not in testfiles_db_filtered[i]:
                        missing_class[i].setdefault(ab_class, []).append(gene)
                        missing_phenotype[i].setdefault(ab_class, {})
                        missing_phenotype[i][ab_class].setdefault(phenotype, []).append(gene)
                        missing_group[i].setdefault(ab_class, {})
                        missing_group[i][ab_class].setdefault(phenotype, {})
                        missing_group[i][ab_class][phenotype].setdefault(group, []).append(gene)
                        missing_gene[i].setdefault(ab_class, {})
                        missing_gene[i][ab_class].setdefault(phenotype, {})
                        missing_gene[i][ab_class][phenotype].setdefault(group, {})
                        missing_gene[i][ab_class][phenotype][group].setdefault(gene, []).append(gene)
                        
                    else:
                        if phenotype not in testfiles_db_filtered[i][ab_class]:
                            missing_phenotype[i].setdefault(ab_class, {})
                            missing_phenotype[i][ab_class].setdefault(phenotype, []).append(gene)
                            missing_group[i].setdefault(ab_class, {})
                            missing_group[i][ab_class].setdefault(phenotype, {})
                            missing_group[i][ab_class][phenotype].setdefault(group, []).append(gene)
                            missing_gene[i].setdefault(ab_class, {})
                            missing_gene[i][ab_class].setdefault(phenotype, {})
                            missing_gene[i][ab_class][phenotype].setdefault(group, {})
                            missing_gene[i][ab_class][phenotype][group].setdefault(gene, []).append(gene)
                            
                        else:
                            if group not in testfiles_db_filtered[i][ab_class][phenotype]:
                                missing_group[i].setdefault(ab_class, {})
                                missing_group[i][ab_class].setdefault(phenotype, {})
                                missing_group[i][ab_class][phenotype].setdefault(group, []).append(gene)
                                missing_gene[i].setdefault(ab_class, {})
                                missing_gene[i][ab_class].setdefault(phenotype, {})
                                missing_gene[i][ab_class][phenotype].setdefault(group, {})
                                missing_gene[i][ab_class][phenotype][group].setdefault(gene, []).append(gene)
                                
                            else:
                                gene_list = testfiles_db_filtered[i][ab_class][phenotype][group]    # Also contains information about identity, coverage and depth
                                if gene not in [gene_list[x][0] for x in range(len(gene_list))]:    # Only gene name and accesion
                                    missing_gene[i].setdefault(ab_class, {})
                                    missing_gene[i][ab_class].setdefault(phenotype, {})
                                    missing_gene[i][ab_class][phenotype].setdefault(group, {})
                                    missing_gene[i][ab_class][phenotype][group].setdefault(gene, []).append(gene)

###########################################################################
# COUNT CORRECT, MISSING AND EXTRA HITS, USED TO GENERATE PLOTS
###########################################################################   
# Count Correct, Extra and Missing hits. Don't include hits with Unknown phenotype.
correct_class_count = list()
missing_class_count = list()
extra_class_count = list() 
correct_phenotype_count = list()
missing_phenotype_count = list()
extra_phenotype_count = list()
correct_group_count = list()
missing_group_count = list()
extra_group_count = list() 
correct_gene_count = list()
missing_gene_count = list()
extra_gene_count = list()
for i in range(len(input_files)):
    # Count for Antibiotic classes  
    correct_class_count.append(len(correct_class[i])) 
    missing_class_count.append(len(missing_class[i])) 
    extra_class_count.append(len(extra_class[i]))
    
    # Count for Phenotypes
    correct_phenotype_count.append(0) 
    missing_phenotype_count.append(0) 
    extra_phenotype_count.append(0)
    for ab_class in correct_phenotype[i]:
        for phenotype in correct_phenotype[i][ab_class]:
            if phenotype != 'Unknown':
                correct_phenotype_count[i] += 1
                
    for ab_class in missing_phenotype[i]:
        for phenotype in missing_phenotype[i][ab_class]:
            if phenotype != 'Unknown':
                missing_phenotype_count[i] += 1   
                
    for ab_class in extra_phenotype[i]:
        for phenotype in extra_phenotype[i][ab_class]:
            if phenotype != 'Unknown':
                extra_phenotype_count[i] += 1
                
    # Count for Gene-groups
    correct_group_count.append(0) 
    missing_group_count.append(0) 
    extra_group_count.append(0)
    for ab_class in correct_group[i]:
        for phenotype in correct_group[i][ab_class]:
            if phenotype != 'Unknown':
                correct_group_count[i] += len(correct_group[i][ab_class][phenotype])
                
    for ab_class in missing_group[i]:
        for phenotype in missing_group[i][ab_class]:
            if phenotype != 'Unknown':
                missing_group_count[i] += len(missing_group[i][ab_class][phenotype])
                
    for ab_class in extra_group[i]:
        for phenotype in extra_group[i][ab_class]:
            if phenotype != 'Unknown':
                extra_group_count[i] += len(extra_group[i][ab_class][phenotype])
                    
    # Count for Gene-variants
    correct_gene_count.append(0) 
    missing_gene_count.append(0) 
    extra_gene_count.append(0)
    for ab_class in correct_gene[i]:
        for phenotype in correct_gene[i][ab_class]:
            if phenotype != 'Unknown':
                for group in correct_gene[i][ab_class][phenotype]:
                    correct_gene_count[i] += len(correct_gene[i][ab_class][phenotype][group])
                    
    for ab_class in missing_gene[i]:
        for phenotype in missing_gene[i][ab_class]:
            if phenotype != 'Unknown':
                for group in missing_gene[i][ab_class][phenotype]:
                    missing_gene_count[i] += len(missing_gene[i][ab_class][phenotype][group])
                    
    for ab_class in extra_gene[i]:
        for phenotype in extra_gene[i][ab_class]:
            if phenotype != 'Unknown':
                for group in extra_gene[i][ab_class][phenotype]:
                    extra_gene_count[i] += len(extra_gene[i][ab_class][phenotype][group])                    

###########################################################################
# OUTPUT RESULTS
########################################################################### 
# Output Counts, used for generating plots
filename = '{}{}_summary_stats.txt'.format(output_path,isolate)
with open(filename, 'w') as outfile:
    outfile.write("Identity cutoff = {}\nCoverage cutoff = {}\nDepth cutoff (Gene-group) = {}%\nMinimum Depth = {}%"
                  .format(args.identity_cutoff,args.coverage_cutoff,100*args.depth_filtration,100*args.min_depth))
    outfile.write('\n##Label\tCorrect hits\tMissing hits\tExtra hits')
    outfile.write('\n\n#Antibiotic class')
    for i in range(len(labels)):
        outfile.write('\n{}\t{}\t{}\t{}'.format(labels[i],correct_class_count[i],missing_class_count[i],extra_class_count[i]))    
    outfile.write('\n\n#Phenotype')
    for i in range(len(labels)):
        outfile.write('\n{}\t{}\t{}\t{}'.format(labels[i],correct_phenotype_count[i],missing_phenotype_count[i],extra_phenotype_count[i])) 
    outfile.write('\n\n#Gene-group')
    for i in range(len(labels)):
        outfile.write('\n{}\t{}\t{}\t{}'.format(labels[i],correct_group_count[i],missing_group_count[i],extra_group_count[i])) 
    outfile.write('\n\n#Gene')
    for i in range(len(labels)):
        outfile.write('\n{}\t{}\t{}\t{}'.format(labels[i],correct_gene_count[i],missing_gene_count[i],extra_gene_count[i])) 


# Output effect of depth filtration
os.chdir(args.input_folder)
input_files = glob.glob('*.res')
input_files.sort()
input_files.pop(-1) # No filtration on Illumina data
for i, file in enumerate(input_files):
    filename = '{}{}_Depth={}_Min_Depth={}_filtration.txt'.format(output_path,file,args.depth_filtration,args.min_depth)
    with open(filename, 'w') as outfile:
        # Output filtration parameters
        outfile.write('Effect of depth filtration on: {}\nCutoffs for including genes are:\nIdentity >= {}\n'
                      'Coverage >= {}\nDepth(Relative to highest depth gene from that gene-group) >= {}%'
                      '\nMinimum gene-group depth (Relative to avg. genome coverage) >= {}%\n'
                      .format(file,args.identity_cutoff,args.coverage_cutoff,100*args.depth_filtration,100*args.min_depth))
        
        # Output results of filtration
        for group in testfiles_gene_groups[i]:
            outfile.write('\nTotal depth of genes from gene-group {} is: {}'.format(group,round(total_depths[i][group],2)))
            # Discarded gene-groups due to min_depth
            if total_depths[i][group] < min_depths[i]:
                outfile.write('\nGene-group not included, due to threshold for minimum gene-group depth\n')
                outfile.write('\nDiscarded {} genes:\n'.format(group))
                for gene in testfiles_gene_groups[i][group]:
                    outfile.write('{} - With a depth of {}\n'.format(gene[0],gene[1][4]))   
            # Included gene-groups
            else:
                max_depth_gene = testfiles_gene_groups[i][group][0]
                outfile.write('\nHighest depth gene from gene-group {} is:\n{} - With a depth of: {}'
                              .format(group,max_depth_gene[0],max_depth_gene[1][4]))
                gene_depth_cutoff = args.depth_filtration*max_depth_gene[1][4]            
                outfile.write('\nDepth threshold for including genes from gene-group {} is: {}'
                              .format(group,round(gene_depth_cutoff,2)))
                # Genes included by depth filtration
                outfile.write('\n\nIncluded {}-genes:\n'.format(group))
                for gene in testfiles_gene_groups[i][group]:
                    gene_depth = gene[1][4]
                    if gene_depth >= gene_depth_cutoff:
                        outfile.write('{} - With a depth of {}\n'.format(gene[0],gene[1][4])) 
                    else:
                        break
                # Genes discarded by depth filtration
                outfile.write('\nDiscarded {}-genes:\n'.format(group))
                for gene in testfiles_gene_groups[i][group]:  
                    gene_depth = gene[1][4]
                    if gene_depth < gene_depth_cutoff:
                        outfile.write('{} - With a depth of {}\n'.format(gene[0],gene[1][4]))
                        
# Output Correct, Missing and Extra hits after filtration, at gene-group level.
os.chdir(args.input_folder)
input_files = glob.glob('*.res')
input_files.sort()
for i, file in enumerate(input_files):
    filename = '{}{}_Depth={}_Min_Depth={}_results.txt'.format(output_path,file,args.depth_filtration,args.min_depth)
    with open(filename, 'w') as outfile:
        outfile.write('Results after depth filtration on: {}\nCutoffs for including genes are:\nIdentity >= {}\n'
                      'Coverage >= {}\nDepth(Relative to highest depth gene from that gene-group) >= {}%'
                      '\nMinimum gene-group depth (Relative to avg. genome coverage) >= {}%\n'
                      .format(file,args.identity_cutoff,args.coverage_cutoff,100*args.depth_filtration,100*args.min_depth))
        outfile.write('Correct, Extra and Missing hits based on results at gene-group level\n')
        outfile.write('Gene information: Gene-variant\tTemplate-Identity\tTemplate-Coverage\tQuery-Identity\tQuery-Coverage\tDepth\n')
        outfile.write('\n############\nCorrect hits\n############')
        for ab_class in correct_group[i]:
            outfile.write('\n###Antimicrobial class - {}'.format(ab_class))
            for phenotype in correct_group[i][ab_class]:
                outfile.write('\n##Phenotype - {}'.format(phenotype))
                for group in correct_group[i][ab_class][phenotype]:
                    genes = correct_group[i][ab_class][phenotype][group]
                    outfile.write('\n#Gene-group - {} - Total depth before filtration, across all phenotypes: {}\n'.format(group,round(total_depths[i][group],2)))
                    outfile.write('Included genes:\n')
                    for gene in genes:
                        gene_name = gene[0]
                        template_identity = gene[1][0]
                        template_coverage = gene[1][1]
                        query_identity = gene[1][2]
                        query_coverage = gene[1][3]
                        depth = gene[1][4]
                        outfile.write('{} - {}\t{}\t{}\t{}\t{}\n'.format(gene_name,template_identity,template_coverage,query_identity,query_coverage,depth))
        outfile.write('\n############\nExtra hits\n############')
        for ab_class in extra_group[i]:
            outfile.write('\n###Antimicrobial class - {}'.format(ab_class))
            for phenotype in extra_group[i][ab_class]:
                outfile.write('\n##Phenotype - {}'.format(phenotype))
                for group in extra_group[i][ab_class][phenotype]:
                    genes = extra_group[i][ab_class][phenotype][group]
                    outfile.write('\n#Gene-group - {} - Total depth before filtration, across all phenotypes: {}\n'.format(group,round(total_depths[i][group],2)))
                    outfile.write('Included genes:\n')
                    for gene in genes:
                        gene_name = gene[0]
                        template_identity = gene[1][0]
                        template_coverage = gene[1][1]
                        query_identity = gene[1][2]
                        query_coverage = gene[1][3]
                        depth = gene[1][4]
                        outfile.write('{} - {}\t{}\t{}\t{}\t{}\n'.format(gene_name,template_identity,template_coverage,query_identity,query_coverage,depth))
        outfile.write('\n############\nMissing hits\n############')
        for ab_class in missing_group[i]:
            outfile.write('\n###Antimicrobial class - {}'.format(ab_class))
            for phenotype in missing_group[i][ab_class]:
                outfile.write('\n##Phenotype - {}\n'.format(phenotype))
                for group in missing_group[i][ab_class][phenotype]:
                    genes = missing_group[i][ab_class][phenotype][group]
                    for gene in genes:
                        gene_name = gene
                        outfile.write('{}\n'.format(gene_name))