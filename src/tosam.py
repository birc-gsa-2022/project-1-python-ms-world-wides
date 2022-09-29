import argparse

def main():
    argparser = argparse.ArgumentParser(description="To Simple-SAM converter")
    argparser.add_argument(
        "mas",
        type=argparse.FileType('r')
    )
    args = argparser.parse_args()

    for line in args.mas:
        chrom, read_name, read_str, pos = line.split('\t')
        
        # Output as Simple-SAM
        print(read_name.strip(), chrom.strip(), str(int(pos)+1), f"{len(read_str.strip())}M", read_str.strip(), sep = '\t')

if __name__ == '__main__':
    main()
