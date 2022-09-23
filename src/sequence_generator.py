import gzip
import random

def cleanup_fasta(fasta_file):
    ''' Function that takes downloaded, gzipped genome fasta file and cleans it
    up to only include the sequence. It then places it into a new text file. This
    can now be used to sample realistic DNA sequences from.'''

    f = gzip.open(fasta_file, 'rt')
    fastafile = f.readlines()

    sequence = []
    i = 0
    for line in fastafile:
        if i == 1000: # change to be able to include more or less
            break
        elif line.startswith('>'):
            continue
        else:
            sequence.append(line.strip().lower())
            i += 1

    final = ''.join(sequence)

    name = 'src/sample_sequence.gz'
    with gzip.open(name, 'wb') as file:
        file.write(final.encode("ascii"))

    return name


def genome_sequence_generator(genome_file, seq_length, s):
    '''Function that takes a gzipped genome file to generate a sequence of length n from
    the genome sequence between the random (with seed s) index i and i+n.'''
    
    f = gzip.open(genome_file, 'rt')
    lines = f.readlines()[0]

    file_len = len(lines)

    if file_len <= seq_length:
        print('Requested sequence length too long')
        return ''
    
    else:
        random.seed(s)

        start = random.randint(0,file_len-seq_length)
        sequence = lines[start:start+seq_length]

        return sequence

    
def random_sequence_generator(length, s, uniform = False, uni_index = 0):
    '''Function that randomize , using seed s, a sequence of bases with given length.
    If uniform is set to True, the sequence will be of a single base according to uni_index.'''
    
    alpha = ['a','t','c', 'g']

    if uniform == True:
        return alpha[uni_index]*length
    else:
        random.seed(s)

        stringlist = [alpha[random.randint(0,3)] for i in range(length)]
        string = ''.join(stringlist)

        return string

def main():
    newname = cleanup_fasta('src/drosophila_melanogaster_genome.gz')

    #newname = 'src/sample_sequence.gz'

    # gsg = genome_sequence_generator(newname, 20, 1)
    # print(gsg)
    # rsg = random_sequence_generator(20, 1, uniform = True, uni_index = 1)
    # print(rsg)


if __name__ == '__main__':
    main()