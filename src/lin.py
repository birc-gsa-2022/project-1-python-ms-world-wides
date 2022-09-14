"""Implementation of a linear time exact matching algorithm."""

import argparse
from fasta_dict import fasta_func
from fastq_dict import fastq_func


def border_algo(x,p):
    array = '$'.join((p,x))
    #built border array
    b = [0]*len(x)


    #Noch zu implementieren!############################
    matches = []

    
    return(matches)

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
    
    #loop through all entries of the dicts and call the naive algorithm
    for p in fastq_dict:
        for x in fasta_dict:
            matches = border_algo(fasta_dict[x], fastq_dict[p])
            for i in matches:
                output(x, p, i, fastq_dict[p])

if __name__ == '__main__':
    main()
