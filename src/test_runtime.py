import time
#import matplotlib.pyplot as plt
from sequence_generator import genome_sequence_generator, random_sequence_generator
from fasta_dict import fasta_func
from fastq_dict import fastq_func
from lin import lin_runner
from naive import naive_runner

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

def change_fasta_length_random(fastqlen, iterations, uniform, matching = True):
    runtime_lin = {}
    runtime_naive = {}

    if matching:
        index = 0
    else:
        index = 1

    fastq_seq = random_sequence_generator(fastqlen, 0, 'fastq', f'read{2}', uniform, uni_index = index)
    # sequence_lengths = range(1000,50000,100)
    sequence_lengths = range(10,50000,1000)

    for j in sequence_lengths:
        fasta_list = []

        for i in range(1,iterations+1):
            l = random_sequence_generator(j, 0, 'fasta', f'seq{i}', uniform) # [name, seq]
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

def plot_runtime(fastqlen1, fastqlen2, iterations, type):
    #plot increasing fasta length

    if type == "compare":

        runtime_lin, runtime_naive = change_fasta_length_random(fastqlen1, iterations, False)
        plt.plot(runtime_lin.keys(), runtime_lin.values(), color="darkblue", label = f"Linear Algorithm with m = {fastqlen1}")
        plt.plot(runtime_naive.keys(), runtime_naive.values(), color="red", label = f"Naive Algorithm with m = {fastqlen1}")

        runtime_lin, runtime_naive = change_fasta_length_random(fastqlen2, iterations, False)
        plt.plot(runtime_lin.keys(), runtime_lin.values(), color="lightblue", label = f"Linear Algorithm with m = {fastqlen2}")
        plt.plot(runtime_naive.keys(), runtime_naive.values(), color="pink", label = f"Naive Algorithm with m = {fastqlen2}")
        plt.xlabel('Size of sequence x [bp]')
        plt.ylabel('Time [s]')
        plt.title(f'Comparing the runtimes of the naive and linear algorithms')
        plt.legend(loc="upper left")
        plt.show()
    else:
        # runtime_lin_rand, runtime_naive_rand = change_fasta_length(fastqlen1, iterations)
        runtime_lin_wc, runtime_naive_wc = change_fasta_length_random(fastqlen1, iterations, True, matching = True)
        runtime_lin_bc, runtime_naive_bc = change_fasta_length_random(fastqlen1, iterations, True, matching = False)
        runtime_lin_rand, runtime_naive_rand = change_fasta_length_random(fastqlen1, iterations, False, matching = False)

        plt.plot(runtime_naive_bc.keys(), runtime_naive_bc.values(),color="#0b8a2b", label = f"Naive Algorithm with m = {fastqlen1}, best case, x = a*n, p = b*m")
        plt.plot(runtime_naive_rand.keys(), runtime_naive_rand.values(), color="#36169e", label = f"Naive Algorithm with m = {fastqlen1} and random sequences")
        plt.plot(runtime_naive_wc.keys(), runtime_naive_wc.values(), color="#d63c75", label = f"Naive Algorithm with m = {fastqlen1}, worst case, x = a*n, p = a*m")
        plt.title(f'Runtime of the naive algorithm')
        plt.xlabel('Size of sequence x [bp]')
        plt.ylabel('Time [s]')
        
        plt.legend(loc="upper left")
        plt.show()
        
        plt.plot(runtime_lin_bc.keys(), runtime_lin_bc.values(),color="#0b8a2b", label = f"Linear Algorithm with m = {fastqlen1}, best case, x = a*n, p = b*m")
        plt.plot(runtime_lin_rand.keys(), runtime_lin_rand.values(), color="#36169e", label = f"Linear Algorithm with m = {fastqlen1} and random sequences")
        plt.plot(runtime_lin_wc.keys(), runtime_lin_wc.values(), color="#d63c75", label = f"Linear Algorithm with m = {fastqlen1}, worst case, x = a*n, p = a*m")
        plt.title(f'Runtime of the linear algorithm')
        plt.xlabel('Size of sequence x [bp]')
        plt.ylabel('Time [s]')

        plt.legend(loc="upper left")
        plt.show()
    

def main():
    plot_runtime(10, 1000, 10, "compare")
    # plot_runtime(10,10, 10, "")

if __name__ == '__main__':
    main()

