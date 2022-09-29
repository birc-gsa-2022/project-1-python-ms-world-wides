import time
from turtle import color
import matplotlib.pyplot as plt
from sequence_generator import genome_sequence_generator, random_sequence_generator
from fasta_dict import fasta_func
from fastq_dict import fastq_func
from lin import lin_runner
from naive import naive_runner
import numpy as np
from math import log

def ntimesm(n,m):
    executionTime_naive = []
    for i in range(n):
        startTime_lin = time.time()
        time.sleep(i*m)
        executionTime_naive.append((time.time() - startTime_lin))
    return executionTime_naive

def nplusm(n,m):
    executionTime_lin = []
    for i in range(n):
        startTime_lin = time.time()
        time.sleep(i+m)
        executionTime_lin.append((time.time() - startTime_lin))
    return executionTime_lin

def change_fasta_length(fastqlen, iterations):
    runtime_lin = {}
    runtime_naive = {}

    fastq_seq = genome_sequence_generator('src/sample_sequence.gz', fastqlen, 0, 'fastq', f'read{2}')
    sequence_lengths = range(10,50000,1000)

    for j in sequence_lengths:
        fasta_list = []

        for i in range(1,iterations+1):
            l = genome_sequence_generator('src/sample_sequence.gz', j, i+j, 'fasta', f'seq{i}') # [name, seq]
            fasta_list += l

        # call function once with multiple fasta
        fasta_d = fasta_func(fasta_list) # dict of multiple fasta, same length
        fastq_d = fastq_func(fastq_seq) # dict of one fastq, set length

        # measure runtime linear algorithm
        startTime_lin = time.time()
        lin_runner(fasta_d, fastq_d)
        executionTime_lin = (time.time() - startTime_lin)
        runtime_lin[j] = executionTime_lin/iterations # add time to value in sequence length key

        # measure runtime naive algorithm
        startTime_naive = time.time()
        naive_runner(fasta_d, fastq_d)
        executionTime_naive = (time.time() - startTime_naive)
        runtime_naive[j] = executionTime_naive/iterations

    return runtime_lin, runtime_naive

def change_fasta_length_uniform(fastqlen, iterations, matching = True):
    runtime_lin = {}
    runtime_naive = {}

    if matching:
        index = 0
    else:
        index = 1

    fastq_seq = random_sequence_generator(fastqlen, 0, 'fastq', f'read{2}', uniform = True, uni_index = index)
    # sequence_lengths = range(1000,50000,100)
    sequence_lengths = range(1000,20000,50)

    for j in sequence_lengths:
        fasta_list = []

        for i in range(1,iterations+1):
            l = random_sequence_generator(j, 0, 'fasta', f'seq{i}', uniform = True) # [name, seq]
            fasta_list += l
        # call function once with multiple fasta
        fasta_d = fasta_func(fasta_list) # dict of multiple fasta, same length
        fastq_d = fastq_func(fastq_seq) # dict of one fastq, set length

        # measure runtime linear algorithm
        startTime_lin = time.time()
        lin_runner(fasta_d, fastq_d)
        executionTime_lin = (time.time() - startTime_lin)
        runtime_lin[j] = executionTime_lin/iterations # add time to value in sequence length key

        # measure runtime naive algorithm
        startTime_naive = time.time()
        naive_runner(fasta_d, fastq_d)
        executionTime_naive = (time.time() - startTime_naive)
        runtime_naive[j] = executionTime_naive/iterations
    
    return runtime_lin, runtime_naive

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

    # ADD PLOT TITLES AND AXIS TITLES

    if type == "compare":

        runtime_lin, runtime_naive = change_fasta_length(fastqlen1, iterations)
        plt.plot(runtime_lin.keys(), runtime_lin.values(), color="darkblue", label = f"Linear Algorithm with m = {fastqlen1}")
        plt.plot(runtime_naive.keys(), runtime_naive.values(), color="red", label = f"Naive Algorithm with m = {fastqlen1}")

        runtime_lin, runtime_naive = change_fasta_length(fastqlen2, iterations)
        plt.plot(runtime_lin.keys(), runtime_lin.values(), color="lightblue", label = f"Linear Algorithm with m = {fastqlen2}")
        plt.plot(runtime_naive.keys(), runtime_naive.values(), color="pink", label = f"Naive Algorithm with m = {fastqlen2}")
    else:
        # runtime_lin_rand, runtime_naive_rand = change_fasta_length(fastqlen1, iterations)
        runtime_lin_wc, runtime_naive_wc = change_fasta_length_uniform(fastqlen1, iterations, matching = False)
        plt.plot(runtime_lin_wc.keys(), runtime_lin_wc.values(), color="darkblue", label = f"Linear Algorithm with m = {fastqlen1}")
        plt.plot(runtime_naive_wc.keys(), runtime_naive_wc.values(), color="red", label = f"Naive Algorithm with m = {fastqlen1}")
        plt.title(f'Runtime of naive and linear algorithms with uniform but not matching sequences')


        # if type == "naive":
        #     plt.plot(x, runtime_naive_bc,color="darkblue", label = "Niave Algorithm with m = 10, x = a^n, p = b^m")
        #     plt.plot(x, runtime_naive_rand, color="red", label = "Naive Algorithm with m = 10 and random sequences")
        #     plt.plot(x, runtime_naive_wc, color="pink", label = "Naive Algorithm with m = 10, x = a^n, p = a^m")

        # if type == "lin":
        #     plt.plot(x, runtime_lin_bc,color="darkblue", label = "Linear Algorithm with m = 10, x = a^n, p = b^m")
        #     plt.plot(x, runtime_lin_rand, color="red", label = "Linear Algorithm with m = 10 and random sequences")
        #     plt.plot(x, runtime_lin_wc, color="pink", label = "Linear Algorithm with m = 10, x = a^n, p = a^m")

    plt.legend(loc="upper left")
    plt.show()
    

def main():
    # plot_runtime(3, 100, 20, "compare")
    plot_runtime(2,10,2, "naive")
    # plot_runtime(10,10,10, "lin")
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
