import time
from matplotlib import pyplot as plt
#
from sequence_generator import genome_sequence_generator
from fasta_dict import fasta_func
from fastq_dict import fastq_func
from lin import lin_runner
from naive import naive_runner

def change_fasta_length(fastqlen, iterations):
    fasta_seq = []
    fastq_seq = []
    runtime_lin = []
    runtime_naive = []


    fastq_seq += genome_sequence_generator('src/sample_sequence.gz', fastqlen, i+5, 'fastq', f'read{i+1}')

    for j in range(0,1000,100):
        for i in range(iterations):
            fasta_seq += genome_sequence_generator('src/sample_sequence.gz', j, i, 'fasta', f'seq{i+1}')

        fasta_d = fasta_func(fasta_seq)
        fastq_d = fastq_func(fastq_seq)
        
        print('Lin:')
        # measure runtime linear algorithm
        startTime_lin = time.time()
        lin_runner(fasta_d, fastq_d)
        executionTime_lin = (time.time() - startTime_lin)
        runtime_lin.append(executionTime_lin)

        print('Naive')
        # measure runtime naive algorithm
        startTime_naive = time.time()
        naive_runner(fasta_d, fastq_d)
        executionTime_naive = (time.time() - startTime_naive)
        runtime_naive.append(executionTime_naive)

    return(range(0,1000,100), runtime_lin, runtime_naive)

def change_fastq_length(fastalen, iterations):
    fasta_seq = []
    fastq_seq = []
    runtime_lin = []
    runtime_naive = []

    fasta_seq += genome_sequence_generator('src/sample_sequence.gz', fastalen, i, 'fasta', f'seq{i+1}')

    for j in range(0,1000,100):
        for i in range(iterations):
            fastq_seq += genome_sequence_generator('src/sample_sequence.gz', j, i+5, 'fastq', f'read{i+1}')

        fasta_d = fasta_func(fasta_seq)
        fastq_d = fastq_func(fastq_seq)
        
        print('Lin:')
        # measure runtime linear algorithm
        startTime_lin = time.time()
        lin_runner(fasta_d, fastq_d)
        executionTime_lin = (time.time() - startTime_lin)
        runtime_lin.append(executionTime_lin)

        print('Naive')
        # measure runtime naive algorithm
        startTime_naive = time.time()
        naive_runner(fasta_d, fastq_d)
        executionTime_naive = (time.time() - startTime_naive)
        runtime_naive.append(executionTime_naive)

    return(range(0,1000,100), runtime_lin, runtime_naive)

def plot_runtime(fastalen, fastqlen, iterations):
    #plot increasing fasta length
    x, runtime_lin, runtime_naive = change_fasta_length(fastqlen, iterations)
    plt.plot(x, runtime_lin)
    plt.plot(x, runtime_naive)
    plt.show()
    
    #plot increasing fastq length
    x, runtime_lin, runtime_naive = change_fastq_length(fastalen, iterations)
    plt.plot(x, runtime_lin)
    plt.plot(x, runtime_naive)
    plt.show()



def main():
    plot_runtime(100, 10, 10)
# test for naive uses no more time than O(nm)


    # check runtime for increasing fasta length



if __name__ == '__main__':
    main()

#check 
# best-case input:  x = A^n y = B^m
# worst-case input: x = A^n y = A^m


# test for lin use no more time than O(n+m)
# Remember to explain your choice of test data.
# What are “best” and “worst” case inputs?
