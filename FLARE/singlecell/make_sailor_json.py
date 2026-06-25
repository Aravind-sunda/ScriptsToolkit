import argparse
import os
from glob import glob
import sys    
import pandas as pd
import json
from pathlib import Path
# Validate paths provided
def validate_paths(paths, verbose=False):
    for to_validate in paths:
        if os.path.exists(to_validate):
            if verbose:
                sys.stdout.write('Confirming {} exists...\n'.format(to_validate))
            else:
                pass
        else:
            sys.stdout.write('Error: {} does not exist!\n'.format(to_validate))
            sys.exit(1)
def create_jsons(samples_path,exp_name, json_folder, sailor_output_folder, fasta, known_snps, edit_type, strand):        
    validate_paths([samples_path, json_folder, sailor_output_folder, fasta], verbose=True)
    samples = glob('{}/*.bam'.format(samples_path))
    #head, tail = os.path.split(samples1)
    #samples = tail
    sys.stdout.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    sys.stdout.write("Found {} finished samples in {}...\n".format(len(samples), samples_path))
    samples1=[]
    for s in samples:
        head, tail= os.path.split(s)
        samples1.append(tail)
        sys.stdout.write('\t{}\n'.format(s))
    sys.stdout.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    sys.stdout.write('Writing JSON inputs for peakcaller for each Samples...\n')
    
    input_json_dicts = []
    
    

    input_json_dict = {
            "samples_path": samples_path,
            "samples": samples1,
            "reverse_stranded": strand,
            "reference_fasta": fasta,
            "known_snps": known_snps,
            "edit_type": edit_type,
            "output_dir": sailor_output_folder
    }

    json_filename = '{}/{}_sailor_input.json'.format(json_folder, exp_name)

    sys.stdout.write("\tWriting {}...\n".format(json_filename))
    with open(json_filename, 'w') as f:
        f.write(json.dumps(input_json_dict, indent=4))

    sys.stdout.write('Done!\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Given finished Bam Index folder, generate appropriate sailor input configuration jsons')
    parser.add_argument('samples_path', type=str)
    parser.add_argument('exp_name', type=str)
    parser.add_argument('json_folder', type=str)
    parser.add_argument('sailor_output_folder', type=str)
    parser.add_argument('fasta', type=str)
    parser.add_argument('known_snps', type=str)
    parser.add_argument('edit_type', type=str)
    parser.add_argument('strand', type=str)

    args = parser.parse_args()
    samples_path = args.samples_path
    exp_name = args.exp_name
    json_folder = args.json_folder
    fasta = args.fasta
    sailor_output_folder = args.sailor_output_folder
    edit_type = args.edit_type
    known_snps = args.known_snps
    strand=args.strand
    create_jsons(samples_path,exp_name, json_folder, sailor_output_folder, fasta, known_snps, edit_type, strand)
