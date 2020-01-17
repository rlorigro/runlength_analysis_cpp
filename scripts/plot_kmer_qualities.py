from collections import defaultdict
from matplotlib import pyplot
import argparse
import sys


def plot_kmer_qualities(path):
    kmer_qualities = list()

    with open(path, "r") as file:
        for l,line in enumerate(file):
            true_kmer, proportion, n = line.split(",")

            proportion = float(proportion)
            n = int(n)

            kmer_qualities.append(proportion)

    fig = pyplot.figure()
    axes = pyplot.axes()
    pyplot.plot(sorted(kmer_qualities))
    axes.set_ylabel("Proportion of matches to reference")
    axes.set_xlabel("Kmer")
    pyplot.show()
    pyplot.close()


def plot_kmer_qualities_vs_n(path):
    kmer_qualities = list()
    n_observations = list()

    with open(path, "r") as file:
        for l,line in enumerate(file):
            true_kmer, proportion, n = line.split(",")

            proportion = float(proportion)
            n = int(n)

            kmer_qualities.append(proportion)
            n_observations.append(n)

    fig = pyplot.figure()
    axes = pyplot.axes()
    pyplot.scatter(n_observations, kmer_qualities, s=0.7, alpha=0.3)
    axes.set_ylabel("Proportion of matches to reference")
    axes.set_xlabel("n")
    pyplot.show()
    pyplot.close()


def main(path):
    plot_kmer_qualities(path)
    plot_kmer_qualities_vs_n(path)


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
