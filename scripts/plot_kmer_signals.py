from matplotlib import pyplot
import os


def read_kmer_signal_means():
    data_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_path = os.path.join(data_path, "data/r9.4_180mv_450bps_6mer.txt")

    means = dict()

    with open(data_path, "r") as file:
        for l,line in enumerate(file):
            if l == 0:
                continue

            data = line.strip().split("\t")
            kmer, level_mean, level_stdv, sd_mean, sd_stdv, weight = data
            level_mean = float(level_mean)

            means[kmer] = level_mean

    return means


def main():
    means = read_kmer_signal_means()

    y = sorted(means.values())
    x = list(range(len(y)))

    pyplot.scatter(x=x,y=y, s=0.3)
    pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    main()