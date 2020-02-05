from matplotlib import pyplot
import argparse
import os


def read_kmer_identities(path):
    identities = dict()

    with open(path, "r") as file:
        for l,line in enumerate(file):
            if l == 0:
                continue

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

            identities[kmer] = n_match/(n_insert + n_delete + n_match + n_mismatch)

    return identities


def read_kmer_signal():
    data_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_path = os.path.join(data_path, "data/r9.4_180mv_450bps_6mer.txt")

    means = dict()
    stdevs = dict()

    with open(data_path, "r") as file:
        for l,line in enumerate(file):
            if l == 0:
                continue

            data = line.strip().split("\t")
            kmer, level_mean, level_stdv, sd_mean, sd_stdv, weight = data
            level_mean = float(level_mean)
            level_stdv = float(level_stdv)

            means[kmer] = level_mean
            stdevs[kmer] = level_stdv

    return means, stdevs


def main(path):
    identities = read_kmer_identities(path)
    means, stdevs = read_kmer_signal()

    axis = pyplot.axes()
    x = list()
    y = list()
    for kmer in identities.keys():
        mean = means[kmer]
        identity = identities[kmer]

        x.append(mean)
        y.append(identity)

    axis.scatter(x=x,y=y, s=0.5, alpha=0.5)
    axis.set_xlabel("Signal mean (pA)")
    axis.set_ylabel("Kmer identity")
    pyplot.show()
    pyplot.close()

    axis = pyplot.axes()
    x = list(range(len(identities)))
    y = list(sorted(identities.values()))

    axis.scatter(x=x,y=y, s=0.3)
    axis.set_xlabel("Kmer")
    axis.set_ylabel("Identity")

    pyplot.show()
    pyplot.close()


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