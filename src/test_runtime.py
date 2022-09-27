import time
from turtle import color
import matplotlib.pyplot as plt
# hello hello 
from sequence_generator import genome_sequence_generator, random_sequence_generator
from fasta_dict import fasta_func
from fastq_dict import fastq_func
from lin import lin_runner
from naive import naive_runner
import numpy as np

def change_fasta_length(fastqlen, iterations):
    fasta_seq = []
    fastq_seq = []
    runtime_lin = []
    runtime_naive = []


    fastq_seq += genome_sequence_generator('src/sample_sequence.gz', fastqlen, 1, 'fastq', f'read{2}')

    for j in range(0,50000,2000):
        for i in range(iterations):
            fasta_seq += genome_sequence_generator('src/sample_sequence.gz', j, i, 'fasta', f'seq{i+1}')

        fasta_d = fasta_func(fasta_seq)
        fastq_d = fastq_func(fastq_seq)

        # measure runtime linear algorithm
        startTime_lin = time.time()
        lin_runner(fasta_d, fastq_d)
        executionTime_lin = (time.time() - startTime_lin)
        runtime_lin.append(executionTime_lin)

        # measure runtime naive algorithm
        startTime_naive = time.time()
        naive_runner(fasta_d, fastq_d)
        executionTime_naive = (time.time() - startTime_naive)
        runtime_naive.append(executionTime_naive)

    return(range(0,50000,2000), runtime_lin, runtime_naive)

def change_fasta_length_uniform(fastqlen, iterations):
    fasta_seq = []
    fastq_seq = []
    runtime_lin = []
    runtime_naive = []

    #length, s, type, name, uniform = False, uni_index = 0
    fastq_seq += random_sequence_generator(fastqlen, )

    for j in range(0,50000,2000):
        for i in range(iterations):
            fasta_seq += random_sequence_generator('src/sample_sequence.gz', j, i, 'fasta', f'seq{i+1}')

        fasta_d = fasta_func(fasta_seq)
        fastq_d = fastq_func(fastq_seq)

        # measure runtime linear algorithm
        startTime_lin = time.time()
        lin_runner(fasta_d, fastq_d)
        executionTime_lin = (time.time() - startTime_lin)
        runtime_lin.append(executionTime_lin)

        # measure runtime naive algorithm
        startTime_naive = time.time()
        naive_runner(fasta_d, fastq_d)
        executionTime_naive = (time.time() - startTime_naive)
        runtime_naive.append(executionTime_naive)

    return(range(0,50000,2000), runtime_lin, runtime_naive)
'''
def change_fastq_length(fastalen, iterations):
    fasta_seq = []
    fastq_seq = []
    runtime_lin = []
    runtime_naive = []

    fasta_seq += genome_sequence_generator('src/sample_sequence.gz', fastalen, 1, 'fasta', f'seq{2}')

    for j in range(0,100,1):
        for i in range(iterations):
            fastq_seq += genome_sequence_generator('src/sample_sequence.gz', j, i+5, 'fastq', f'read{i+1}')

        fasta_d = fasta_func(fasta_seq)
        fastq_d = fastq_func(fastq_seq)

        # measure runtime linear algorithm
        startTime_lin = time.time()
        lin_runner(fasta_d, fastq_d)
        executionTime_lin = (time.time() - startTime_lin)
        runtime_lin.append(executionTime_lin)

        # measure runtime naive algorithm
        startTime_naive = time.time()
        naive_runner(fasta_d, fastq_d)
        executionTime_naive = (time.time() - startTime_naive)
        runtime_naive.append(executionTime_naive)

    return(range(0,100,1), runtime_lin, runtime_naive)
'''
def plot_runtime(fastqlen1, fastqlen2, iterations, type):
    #plot increasing fasta length

    if type == "compare":
        x, runtime_lin, runtime_naive = change_fasta_length(fastqlen1, iterations)
        plt.plot(x, runtime_lin,color="darkblue", label = "Linear Algorithm with m = 3")
        plt.plot(x, runtime_naive, colorNaive="red", label = "Naive Algorithm with m = 3")

        x, runtime_lin, runtime_naive = change_fasta_length(fastqlen2, iterations)
        plt.plot(x, runtime_lin, color="lightblue", label = "Linear Algorithm with m = 100")
        plt.plot(x, runtime_naive, color="pink", label = "Naive Algorithm with m = 100")
    else:
        x, runtime_lin_rand, runtime_naive_rand = change_fasta_length(fastqlen1, iterations)
        x, runtime_naive_wc, runtime_lin_wc= change_fasta_length_uniform(fastqlen1, iterations)
        if type == "naive":
            plt.plot(x, runtime_naive_bc,color="darkblue", label = "Niave Algorithm with m = 10, x = a^n, p = b^m")
            plt.plot(x, runtime_naive_rand, color="red", label = "Naive Algorithm with m = 10 and random sequences")
            plt.plot(x, runtime_naive_wc, color="pink", label = "Naive Algorithm with m = 10, x = a^n, p = a^m")

        if type == "lin":
            plt.plot(x, runtime_lin_bc,color="darkblue", label = "Linear Algorithm with m = 10, x = a^n, p = b^m")
            plt.plot(x, runtime_lin_rand, color="red", label = "Linear Algorithm with m = 10 and random sequences")
            plt.plot(x, runtime_lin_wc, color="pink", label = "Linear Algorithm with m = 10, x = a^n, p = a^m")

    plt.legend(loc="upper left")
    plt.show()
    

def main():
    plot_runtime(3, 100, 10, "compare")
    plot_runtime(10,10,10, "naive")
    plot_runtime(10,10,10, "lin")
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
