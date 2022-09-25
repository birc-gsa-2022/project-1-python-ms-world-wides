from sequence_generator import genome_sequence_generator
from fasta_dict import fasta_func
from fastq_dict import fastq_func
from lin import lin_runner
from naive import naive_runner
#try
# test for naive uses no more time than O(nm)

# have the same seed as length?

def main():
    fasta_seq = genome_sequence_generator('src/sample_sequence.gz', 50, 1, 'fasta')
    fastq_seq = genome_sequence_generator('src/sample_sequence.gz', 2, 2, 'fastq')

    fasta_d = fasta_func(fasta_seq)
    fastq_d = fastq_func(fastq_seq)
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
