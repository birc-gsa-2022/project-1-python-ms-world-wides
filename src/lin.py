"""Implementation of a linear time exact matching algorithm."""

import argparse
from fasta_dict import fasta_func
from fastq_dict import fastq_func

def border_algo(x,p):
    #edge case
    if p == "" or x == "":
        return []
    #create string
    jointSeq = '$'.join((p,x))
    #build border array
    ba = [0]*len(jointSeq)
    
    p_len = len(p)
    result = []
    #otherwise run trough seq
    for i in range(1, len(jointSeq)):
        b = ba[i-1]
        #extend border
        while (b > 0) and (jointSeq[i] != jointSeq[b]):
            b = ba[b-1]
        #if matches
        if(jointSeq[i]==jointSeq[b]):
            ba[i] = b + 1
            # Add to result already in the first loop to reduce runtime
            if b + 1 == p_len and i >= p_len:
                result.append(i-2*len(p)+1)
        else:
            ba[i] = 0

    return result


#prints the data in a simple sam-format
def output(x_name, p_name, i, p):
    return 

def lin_runner(fasta_dict, fastq_dict):

    for p_key, p_val in fastq_dict.items():
        for x_key, x_val in fasta_dict.items():
            matches = border_algo(x_val, p_val)
            for i in matches:
                print('\t'.join([x_key, p_key, str(i), f'{str(len(p_val))}M', p_val]))

def main():
    
    argparser = argparse.ArgumentParser(
        description="Exact matching in linear time")
    argparser.add_argument("genome", type=argparse.FileType('r'))
    argparser.add_argument("reads", type=argparse.FileType('r'))
    args = argparser.parse_args()

   #translate files into dicts
    fasta_dict = fasta_func(args.genome)
    fastq_dict = fastq_func(args.reads)

    lin_runner(fasta_dict, fastq_dict)

if __name__ == '__main__':
    main()
