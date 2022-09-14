import argparse

def fasta_func(fastafile):

    sequence = []
    name = ''
    fasta_dict = {}
    for line in fastafile:
        if line.startswith('>'):
            if name != '':
                fasta_dict[name] = ''.join(sequence)
                sequence = []
            name = line[1:].strip()
        else:
            sequence.append(line.strip())

    if name != '':
        fasta_dict[name] = ''.join(sequence)

    return fasta_dict


def main():
    argparser = argparse.ArgumentParser(
        description="Extract Simple-FASTA records"
    )
    argparser.add_argument(
        "fasta",
        type=argparse.FileType('r')
    )
    args = argparser.parse_args()
    dic = fasta_func(args.fasta)
    print(dic)

if __name__ == '__main__':
    main()
