from collections import defaultdict,Counter
from matplotlib import pyplot
import argparse

def get_length_distributions(path):
    operations = defaultdict(lambda:defaultdict(int))
    n_string = ""

    with open(path) as file:
        for l,line in enumerate(file):
            # Skip all header lines
            if line.startswith("@"):
                continue

            # Fetch the cigar string
            cigar_string = line.split("\t")[5]

            # Skip mappings with no alignment (0 MQ)
            if cigar_string == "*":
                continue

            # Iterate the cigar string
            for c in cigar_string:
                if c.isnumeric():
                    n_string += c
                else:
                    # Increment the operation in the dictionary
                    operations[c][int(n_string)] += 1
                    n_string = ""

    for item in operations.items():
        if item[0] in {"M","X","I","D","="}:
            frequencies = list(sorted(item[1].items(),key=lambda x: x[0]))
            x,y = zip(*frequencies)

            pyplot.plot(x,y)
            pyplot.title(item[0])
            pyplot.show()
            pyplot.close()


def main(path):
    operations = defaultdict(int)
    n_string = ""

    with open(path) as file:
        for l,line in enumerate(file):
            # Skip all header lines
            if line.startswith("@"):
                continue

            # Fetch the cigar string
            line = line.split("\t")
            cigar_string = line[5]
            flag = int(line[1])

            # only allow F/R primary, supplementary
            if flag not in {0,16,2064,2048}:
                continue

            # Skip mappings with no alignment (0 MQ)
            if cigar_string == "*":
                continue

            # Iterate the cigar string
            for c in cigar_string:
                if c.isnumeric():
                    n_string += c
                else:
                    n = int(n_string)

                    # Skip SVs
                    if c in {"I", "D"}:
                        # Skip large indels (SVs?)
                        if n > 500:
                            n_string = ""
                            continue
                        else:
                            # Increment the operation in the dictionary
                            operations[c] += n
                            n_string = ""
                    else:
                        # Increment the operation in the dictionary
                        operations[c] += n
                        n_string = ""

    print(operations)

    sum = 0
    for item in operations.items():
        if item[0] in {"X","I","D","="}:
            print("%s: %d"%(item[0], item[1]))
            sum += item[1]

    for item in operations.items():
        if item[0] in {"X","I","D","="}:
            print("%s rate: %f"%(item[0], item[1]/sum))

    for item in operations.items():
        if item[0] in {"X","I","D","="}:
            print("%s per 100kb: %f"%(item[0], (item[1]/sum)*100_000))


if __name__ == "__main__":
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="FASTA file path of true reference to be compared against"
    )

    args = parser.parse_args()

    main(path=args.bam)
    get_length_distributions(path=args.bam)
