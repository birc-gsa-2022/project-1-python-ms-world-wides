import argparse
from fasta_dict import fasta_func
from fastq_dict import fastq_func

def naive_align(x, p):

    match_indexes = []
    for i in range(len(x)-len(p)+1):
        for j, el_p in enumerate(p):
            if x[i+j]!=el_p:
                break
        else:
            match_indexes.append(i+1)
    
    return match_indexes

def output(x_name, p_name, i, p):
    '''Function that takes the name of sequence x to which the pattern p (p_name) has been exactly
    aligned at index i. The function prints the data in a simple sam-format'''
    
    return print(p_name, x_name, i, f'{str(len(p))}M', p, sep = '\t')
    


def main():
    # OBS only works with exact matches
    # help
    argparser = argparse.ArgumentParser(
        description="Extract Simple-FASTA and Simple-FASTQrecords"
    )
    argparser.add_argument("fasta", type=argparse.FileType('r'))
    argparser.add_argument("fastq", type=argparse.FileType('r'))

    args = argparser.parse_args()

    fasta_dict = fasta_func(args.fasta)
    fastq_dict = fastq_func(args.fastq)

    for kp, vp in fastq_dict.items():
        for kx, vx in fasta_dict.items():
            index_list = naive_align(vx, vp)
            for i in index_list:
                output(kx, kp, i, vp)

if __name__ == '__main__':
    main()