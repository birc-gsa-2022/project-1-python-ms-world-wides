import argparse
import os
from fasta_dict import fasta_func
from fastq_dict import fastq_func

# input: x from fasta_dict.py
# input: p from fast1_dict.py

# script: align p to x using the naive algorithm. 
# return name(x), name(p), start index (from 1), cigar (M only), p.

# functions: alignment function, returning index
# cigar function called on in alignment function?

# output: simple_sam format. Use pd dataframe maybe?

def naive_align(x, p):

    match_indexes = []
    for i in x:
        for j in p:
            if j == len(p):
                if x[i+j] == p[j]:
                    match_indexes.append(i+1) # indexing from 1
                else:
                    break
            elif x[i+j] == p[j]:
                continue
            else:
                break
    
    return match_indexes

def output(x_name, p_name, i, p):
    '''Function that takes the name of sequence x to which the pattern p (p_name) has been exactly
    aligned at index i. The function prints the data in a simple sam-format'''
    
    return print(p_name, x_name, i, f'{str(len(p))}M', p, sep = '\t')


def main():

    # fasta files here
    argparser = argparse.ArgumentParser(
        description="Extract Simple-FASTA and Simple-FASTQrecords"
    )
    argparser.add_argument("fasta", type=argparse.FileType('r'))
    argparser.add_argument("fastq", type=argparse.FileType('r'))

    args = argparser.parse_args()

    fasta_dict = fasta_func(args.fasta)
    fastq_dict = fastq_func(args.fastq)

    pass


if __name__ == '__main__':
    main()