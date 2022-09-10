import argparse

def main():
    argparser = argparse.ArgumentParser(
        description="Extract Simple-FASTA records"
    )
    argparser.add_argument(
        "fasta",
        type=argparse.FileType('r')
    )
    args = argparser.parse_args()

    # print(f"Now I need to process the records in {args.fasta}")

    printer = []
    for line in args.fasta:
        if line.startswith('>'):
            if len(printer) != 0:
                print(printer)
                printer = ''
            name = line[1:]
            name = name.strip()
            printer.append(name)
            printer.append('\t')
        else:
            printer.append(line.strip())

    printer = ''.join(printer)
    # print last line     
    print(printer)

if __name__ == '__main__':
    main()
