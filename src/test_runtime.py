from re import I
from sequence_generator import genome_sequence_generator
from fasta_dict import fasta_func
from fastq_dict import fastq_func
from lin import lin_runner
from naive import naive_runner
#try
# test for naive uses no more time than O(nm)

# have the same seed as length?

def main():
    fasta_seq = []
    fastq_seq = []
    for i in range(2):
        fasta_seq += genome_sequence_generator('src/sample_sequence.gz', 50, i, 'fasta', f'seq{i+1}')
        fastq_seq += genome_sequence_generator('src/sample_sequence.gz', 3, i+5, 'fastq', f'read{i+1}')

    print(fasta_seq)
    print(fastq_seq)

    fasta_d = fasta_func(fasta_seq)
    print(fasta_d)
    fastq_d = fastq_func(fastq_seq)
    print(fastq_d)
    print('Lin:')
    print(lin_runner(fasta_d, fastq_d))
    print('Naive')
    print(naive_runner(fasta_d, fastq_d))

    return 0

if __name__ == '__main__':
    main()

#check 
# Remember to explain your choice of test data. 


# best-case input:  x = A^n y = B^m
# worst-case input: x = A^n y = A^m


# test for lin use no more time than O(n+m)
# Remember to explain your choice of test data.
# What are “best” and “worst” case inputs?
