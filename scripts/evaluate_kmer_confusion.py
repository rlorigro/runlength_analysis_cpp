from collections import defaultdict
import argparse
import sys


def validate_number_of_true_kmers(path):
    true_kmers = set()

    with open(path, "r") as file:
        for l,line in enumerate(file):
            line = line.split(",")
            true_kmers.add(line[0])

            if l % 1_000_000 == 0:
                print(l, len(true_kmers))

    print(len(true_kmers))


def find_ref_kmer_qualities(path):
    kmer_qualities = defaultdict(list)

    with open(path, "r") as file:
        for l,line in enumerate(file):
            true_kmer, observed_kmer, proportion, n = line.split(",")

            proportion = float(proportion)
            n = int(n)

            # Find the matches, ignore mismatches
            if true_kmer == observed_kmer:
                kmer_qualities[true_kmer] = [proportion, n]

            if l % 1_000_000 == 0:
                sys.stderr.write(str(l) + " " + str(len(kmer_qualities)) + "\n")
                sys.stderr.flush()

    for item in sorted(kmer_qualities.items(), key=lambda x: x[1][0]):
        data = [item[0], item[1][0], item[1][1]]
        print(",".join(list(map(str,data))))


def find_read_kmer_qualities(path):
    kmer_qualities = defaultdict(list)

    with open(path, "r") as file:
        for l,line in enumerate(file):
            true_kmer, observed_kmer, proportion, n_ref = line.split(",")

            proportion = float(proportion)
            n_ref = int(n_ref)

            # Find the matches, ignore mismatches
            if true_kmer == observed_kmer:
                kmer_qualities[true_kmer] = [proportion, n_ref]

            if l % 1_000_000 == 0:
                sys.stderr.write(str(l) + " " + str(len(kmer_qualities)) + "\n")
                sys.stderr.flush()

    for item in sorted(kmer_qualities.items(), key=lambda x: x[1][0]):
        data = [item[0], item[1][0], item[1][1]]
        print(",".join(list(map(str,data))))


def print_kmer_quality_range(path):
    with open(path, "r") as file:
        for l,line in enumerate(file):
            true_kmer, proportion, n = line.split(",")

            proportion = float(proportion)
            n = int(n)

            if 0.600 < proportion < 0.670:
                print(true_kmer)


def main(path):
    # find_kmer_qualities(path)
    print_kmer_quality_range(path)


if __name__ == "__main__":
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="csv of kmer confusion rates"
    )

    args = parser.parse_args()

    main(path=args.input)
