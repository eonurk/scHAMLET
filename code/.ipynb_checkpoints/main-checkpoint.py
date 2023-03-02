import os
import sys
import json
import pysam
import subprocess
from collections import Counter

# activate venv/ environment

def run_freebayes(bam_file, reference_file, freebayes_vcf_output, important_variations_file):

    command = "freebayes -f " + reference_file + " -C 5 " +  bam_file + " -t " + important_variations_file + " > " + freebayes_vcf_output
    print("\n", command)
    process = subprocess.Popen(command, shell=True)
    process.wait()
    print("Variant calling done.")

def index_vcf(vcf_file):
    command = "bgzip -c " + vcf_file + " > " + vcf_file + ".gz"
    process = subprocess.Popen(command, shell=True)
    process.wait()
    command = "tabix -p vcf " + vcf_file +".gz"
    process = subprocess.Popen(command, shell=True)
    print("Indexing vcf file done.")
    process.wait()


def read_important_variations(filename):
    important_variations = []
    with open(filename) as important_variations_file:
        for line in important_variations_file:
            line = line.rstrip()
            line = line.split("\t")
            line[1] = int(line[1])
            line[2] = int(line[2])
            important_variations.append(line[0:3])
    return important_variations

def printVars(object):
    for i in [v for v in dir(object) if not callable(getattr(object,v))]:
        print('\n%s:' % i)
        print(getattr(object, i))

def variations_in_vcf(filepath, important_variations):
    vcf_file = pysam.VariantFile(filepath)
    existing_imp_var = []
    for rec in vcf_file.fetch():
        if important_variations != None:
            # check later if rec.alts is > 0
            existing_imp_var.append([rec.contig, rec.start, rec.stop, rec.ref, rec.alts[0]]) 
        else:
            existing_imp_var.append(rec)
    print(len(existing_imp_var), " variations in total")
    return existing_imp_var

def call_snp(imp_var, alignment_file):
    cell_barcodes = {'wt': set(), 'mut': set()}
    alignment = alignment_file.fetch(imp_var[0], imp_var[1], imp_var[2])
    
    for read in alignment:
        for pair in read.aligned_pairs:
            
            if (pair[0] == None or pair[1] == None):
                continue
            
            if (pair[1] == imp_var[1]): # be careful with indexing vcf vs. calling 
                
                if (read.query_sequence[pair[0]] == imp_var[4]):
                    try:        
                        # if snp get read cb 
                        cell_barcodes['mut'].add(read.get_tag("CB"))     
                    except:
                        continue
                else:
                    try:
                        # if no indel then wild type
                        cell_barcodes['wt'].add(read.get_tag("CB"))
                    except:
                        continue
                    
    return cell_barcodes

def call_indel(imp_var, alignment_file):
    
    cell_barcodes = {'wt': set(), 'mut': set()}
    for pileupcolumn in alignment_file.pileup(imp_var[0], imp_var[1], imp_var[2]):

        # limit search region to the variant 
        if (pileupcolumn.pos<= imp_var[2]) and (pileupcolumn.pos >= imp_var[1]):          

            for pileupread in pileupcolumn.pileups:

                # check indels
                if pileupread.indel != 0:

                    try:        
                        # if indel add mut
                        cell_barcodes['mut'].add(pileupread.alignment.get_tag("CB"))     
                    except:
                        continue
                else:
                    try:
                        # if no indel then wild type
                        cell_barcodes['wt'].add(pileupread.alignment.get_tag("CB"))
                    except:
                        continue
    return cell_barcodes

def find_reads_with_variation(bam_file, existing_imp_var):
    
    
    alignment_file = pysam.AlignmentFile(bam_file, "rb", check_sq=False)
    
    var_dict = {}
    for imp_var in existing_imp_var:
        
        if len(imp_var[4]) == len(imp_var[3]):  # SNPs 
            cell_barcodes = call_snp(imp_var, alignment_file)
        else:
            cell_barcodes = call_indel(imp_var, alignment_file)
            
        print("variation is ", imp_var)
                        
        var_key = imp_var[0]+"_"+str(imp_var[1])+"_"+str(imp_var[2])+"_"+imp_var[3]+"_"+imp_var[4]
        
        # remove the mut cells that are assigned to wt         
        cell_barcodes['wt'] = cell_barcodes['wt'] - cell_barcodes['mut']
        
        cell_barcodes['wt'] = list(cell_barcodes['wt']) 
        cell_barcodes['mut'] = list(cell_barcodes['mut'])
        
        var_dict[var_key] = cell_barcodes
        
    return var_dict
    

#Press the green button in the gutter to run the script.
if __name__ == '__main__':
    
    bam_file_txt = sys.argv[1]
    important_variations_file = sys.argv[2]
    reference_file = sys.argv[3]
    output_vcf_file = ""
    
    # Using readlines()
    bam_locs = open(bam_file_txt, 'r')
    bam_files = bam_locs.readlines()
    
    if not os.path.exists("vcf"):
        os.mkdir("vcf")
        print("vcf directory created")
    
    if important_variations_file == "None":
        print("No variation found")
    
    else:
        
        important_variations = read_important_variations(important_variations_file)
        for bam_file in bam_files:
            
            bam_file = bam_file.strip()
            sample_name = bam_file.split("/")[8].replace("_GEX_count", "").replace("_count", "")
            
            if (sample_name == "GBZ_10195"):
                sample_name = "scAML_1_003_10197"
            elif(sample_name == "LAE_5904"):
                sample_name = "scAML_3_005_5904"
                
            print("\n" + sample_name)
            
            output_vcf_file = "vcf/" + sample_name + ".vcf"
            
            run_freebayes(bam_file, reference_file, output_vcf_file, important_variations_file)
            index_vcf(output_vcf_file)
            existing_imp_var = variations_in_vcf(output_vcf_file + ".gz", important_variations)
            var_dict = find_reads_with_variation(bam_file, existing_imp_var)
            
            with open("output/" + sample_name + ".json", "w") as write_file:
                json.dump(var_dict, write_file)