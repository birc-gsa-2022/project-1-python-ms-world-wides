"""Implementation of a linear time exact matching algorithm."""

import argparse
from ast import pattern
from fasta_dict import fasta_func
from fastq_dict import fastq_func

#
def border_algo(x,p):
    #edge case
    if p == "" or x == "":
        return []
    #crate string
    jointSeq = '$'.join((p,x))
    #built border array
    ba = [0]*len(jointSeq)
    
    #otherwise run trough seq
    for i in range(1, len(jointSeq)):
        b = ba[i-1]
        #extend border
        while (b > 0) and (jointSeq[i] != jointSeq[b]):
            b = ba[b-1]
        #if matches
        if(jointSeq[i]==jointSeq[b]):
            ba[i] = b + 1
        else:
            ba[i] = 0

    #filter matches with length of pattern
    result = []
    for i in range(len(p)-1, len(ba)):
        if ba[i] == len(p):
            result.append(i-len(p)-1)
    return result

#prints the data in a simple sam-format
def output(x_name, p_name, i, p):
    return print(p_name, x_name, i, f'{str(len(p))}M', p, sep = '\t')

def main():
    argparser = argparse.ArgumentParser(
        description="Exact matching in linear time")
    argparser.add_argument("genome", type=argparse.FileType('r'))
    argparser.add_argument("reads", type=argparse.FileType('r'))
    args = argparser.parse_args()

   #translate files into dicts
    fasta_dict = fasta_func(args.genome)
    fastq_dict = fastq_func(args.reads)
    
    #Test
    # print(border_algo('ABABABAA','AB'))
    # print(border_algo('ABABABAA',''))
    # print(border_algo('',''))
    # print(border_algo('AAAA','B'))
    # print(border_algo('AAAA','A'))

    #loop through all entries of the dicts and call the naive algorithm
    for p in fastq_dict:
        for x in fasta_dict:
            matches = border_algo(fasta_dict[x], fastq_dict[p])
            for i in matches:
                output(x, p, i, fastq_dict[p])

if __name__ == '__main__':
    main()
