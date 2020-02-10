from collections import defaultdict
from matplotlib import pyplot
import argparse


class CigarCount:
    def __init__(self):
        self.n_align_match = 0
        self.n_insert = 0
        self.n_delete = 0
        self.n_refskip = 0
        self.n_softclip = 0
        self.n_hardclip = 0
        self.n_pad = 0
        self.n_match = 0
        self.n_mismatch = 0

    def __str__(self):
        s = str(self.n_align_match) + "," + \
            str(self.n_insert) + "," + \
            str(self.n_delete) + "," + \
            str(self.n_refskip) + "," + \
            str(self.n_softclip) + "," + \
            str(self.n_hardclip) + "," + \
            str(self.n_pad) + "," + \
            str(self.n_match) + "," + \
            str(self.n_mismatch) + ","
        return s


def complement_base(base):
    if base == "A":
        return "T"
    elif base == "C":
        return "G"
    elif base == "G":
        return "C"
    elif base == "T":
        return "A"
    else:
        exit("ERROR: invalid base has no complement: " + base)


def reverse_complement(kmer):
    rc = ""
    for base in reversed(kmer):
        rc += complement_base(base)

    return rc


def aggregate_reverse_complement_cigar_counts(path):
    bidirectional_kmer_cigar_counts = defaultdict(lambda: CigarCount())

    with open(path, "r") as file:
        for l,line in enumerate(file):
            data = line.strip().split(",")
            kmer, n_align_match, n_insert, n_delete, n_refskip, n_softclip, n_hardclip, n_pad, n_match, n_mismatch, _ = data

            n_align_match = int(n_align_match)
            n_insert = int(n_insert)
            n_delete = int(n_delete)
            n_refskip = int(n_refskip)
            n_softclip = int(n_softclip)
            n_hardclip = int(n_hardclip)
            n_pad = int(n_pad)
            n_match = int(n_match)
            n_mismatch = int(n_mismatch)

            kmer_complement = reverse_complement(kmer)

            min_kmer = min(kmer, kmer_complement)

            bidirectional_kmer_cigar_counts[min_kmer].n_align_match += n_align_match
            bidirectional_kmer_cigar_counts[min_kmer].n_insert += n_insert
            bidirectional_kmer_cigar_counts[min_kmer].n_delete += n_delete
            bidirectional_kmer_cigar_counts[min_kmer].n_refskip += n_refskip
            bidirectional_kmer_cigar_counts[min_kmer].n_softclip += n_softclip
            bidirectional_kmer_cigar_counts[min_kmer].n_hardclip += n_hardclip
            bidirectional_kmer_cigar_counts[min_kmer].n_pad += n_pad
            bidirectional_kmer_cigar_counts[min_kmer].n_match += n_match
            bidirectional_kmer_cigar_counts[min_kmer].n_mismatch += n_mismatch

            # print("--")
            # print(line.strip())
            # print(kmer, min_kmer, str(bidirectional_kmer_cigar_counts[min_kmer]))

    return bidirectional_kmer_cigar_counts


def main(path):
    output_path = ".".join(path.strip().split(".")[:-1]) + "_bidirectional.csv"
    print("Writing file: " + output_path)

    bidirectional_kmer_cigar_counts = aggregate_reverse_complement_cigar_counts(path)

    with open(output_path, "w") as file:
        for kmer, cigar_counts in bidirectional_kmer_cigar_counts.items():
            file.write(kmer)
            file.write(",")
            file.write(str(cigar_counts))
            file.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="path of file containing quadrant bounds"
    )
    args = parser.parse_args()

    main(path=args.input)