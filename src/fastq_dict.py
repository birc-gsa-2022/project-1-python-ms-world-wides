import argparse

def main():
    argparser = argparse.ArgumentParser(
        description="Extract the sequences from a simple-fastq file")
    argparser.add_argument(
        "fastq",
        type=argparse.FileType('r')
    )
    args = argparser.parse_args()

    read = []
    name = ''
    fastq_dict = {}
    for line in args.fastq:
        if line.startswith('@'):
            if name != '':
                fastq_dict[name] = ''.join(read)
                read = []
            name = line[1:].strip()
        else:
            read.append(line.strip())

    if name != '':
        fastq_dict[name] = ''.join(read)

    return fastq_dict

if __name__ == '__main__':
    main()
