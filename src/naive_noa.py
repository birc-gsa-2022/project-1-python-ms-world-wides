"""Implementation of the naive exact matching algorithm."""

import argparse


def main():
    argparser = argparse.ArgumentParser(
        description="Exact matching using the naive method")
    argparser.add_argument("genome", type=argparse.FileType('r'))
    argparser.add_argument("reads", type=argparse.FileType('r'))
    args = argparser.parse_args()
    print(f"Find every reads in {args.reads.name} " +
          f"in genome {args.genome.name}")


    list = []

    for g in args.genome:

        n = len(g.sequence)

        for r in args.reads:

            m = len(args.reads)
            
            for i in range(n - m + 1):
                j = 0
                while j<m :
                    if (g.sequence[i+j]!=r.read[j]):
                        break
                
                if (j==m):
                    list.append([r.name, g.name, i, m+"M", r.read])

    return list



if __name__ == '__main__':
    main()
