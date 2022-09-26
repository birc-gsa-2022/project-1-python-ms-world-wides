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
    
def naive_runner(fasta_dict, fastq_dict):

    l = []
    for p_key, p_val in fastq_dict.items():
        for x_key, x_val in fasta_dict.items():
            matches = naive_align(x_val, p_val)
            for i in matches:
                l.append('\t'.join([p_key, x_key, str(i), f'{str(len(p_val))}M', p_val]))
    
    return '\n'.join(l)

def main():

    argparser = argparse.ArgumentParser(
        description="Extract Simple-FASTA and Simple-FASTQrecords"
    )
    argparser.add_argument("genome", type=argparse.FileType('r'))
    argparser.add_argument("reads", type=argparse.FileType('r'))

    args = argparser.parse_args()

    fasta_dict = fasta_func(args.genome)
    fastq_dict = fastq_func(args.reads)

    print(naive_runner(fasta_dict, fastq_dict))

if __name__ == '__main__':
    main()