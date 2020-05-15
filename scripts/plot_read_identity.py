import argparse
from matplotlib import pyplot
import numpy


def plot_2d(paths):
    colors = [(0.945, 0.267, 0.176), (0.122, 0.498, 0.584), (0.945, 0.71, 0.176)]

    fig = pyplot.figure()
    axes = pyplot.axes()

    labels = list()

    for p,path in enumerate(paths):
        lengths = list()
        identities = list()

        with open(path, "r") as file:
            for line in file:
                data = line.strip().split(",")
                id, n_match, n_mismatch, n_insert, n_delete = data

                n_match = int(n_match)
                n_mismatch = int(n_mismatch)
                n_insert = int(n_insert)
                n_delete = int(n_delete)

                if n_insert + n_delete + n_match + n_mismatch == 0:
                    continue

                identity = n_match/(n_insert + n_delete + n_match + n_mismatch)
                length = sum([n_insert + n_delete + n_match + n_mismatch])
                lengths.append(length)
                identities.append(identity)

        axes.scatter(x=lengths,y=identities, color=colors[p],s=0.3,alpha=0.3)

        axes.set_xlabel("Length")
        axes.set_ylabel("Proportion")

        name = path.split("/")[-1][26:].split("_VS")[0]
        print(name)
        labels.append(name)

    axes.legend(labels)

    fig.set_size_inches(12,6)
    pyplot.savefig("Read_identities.png",dpi=200)

    pyplot.show()
    pyplot.close()


def summary_stats(paths):
    for p,path in enumerate(paths):
        identities = list()
        mismatch_rates = list()
        insert_rates = list()
        delete_rates = list()

        with open(path, "r") as file:
            for line in file:
                data = line.strip().split(",")
                id, n_match, n_mismatch, n_insert, n_delete = data

                n_match = int(n_match)
                n_mismatch = int(n_mismatch)
                n_insert = int(n_insert)
                n_delete = int(n_delete)

                if n_insert + n_delete + n_match + n_mismatch == 0:
                    continue

                if n_insert + n_delete + n_match + n_mismatch < 10000:
                    continue

                identity = n_match/(n_insert + n_delete + n_match + n_mismatch)
                mismatch_rate = n_mismatch/(n_insert + n_delete + n_match + n_mismatch)
                insert_rate = n_insert/(n_insert + n_delete + n_match + n_mismatch)
                delete_rate = n_delete/(n_insert + n_delete + n_match + n_mismatch)

                identities.append(identity)
                mismatch_rates.append(mismatch_rate)
                insert_rates.append(insert_rate)
                delete_rates.append(delete_rate)

        print(path)
        print("median_identity\t" + str(numpy.median(numpy.array(identities))))
        print("median_mismatch_rate\t" + str(numpy.median(numpy.array(mismatch_rates))))
        print("median_insert_rate\t" + str(numpy.median(numpy.array(insert_rates))))
        print("median_delete_rate\t" + str(numpy.median(numpy.array(delete_rates))))



def main(paths):
    colors = [(0.945, 0.267, 0.176),
              (0.122, 0.498, 0.584),
              (0.945, 0.71, 0.176),
              (0.58,0.761,0.153),
              (0.788,0.161,0.208),
              (0.38,0.137,0.545)]

    fig = pyplot.figure()
    axes = pyplot.axes()

    labels = list()

    for p,path in enumerate(paths):
        identities = list()

        with open(path, "r") as file:
            for line in file:
                data = line.strip().split(",")
                id, n_match, n_mismatch, n_insert, n_delete = data

                n_match = int(n_match)
                n_mismatch = int(n_mismatch)
                n_insert = int(n_insert)
                n_delete = int(n_delete)

                if n_insert + n_delete + n_match + n_mismatch == 0:
                    continue

                if n_insert + n_delete + n_match + n_mismatch < 10000:
                    continue

                identity = n_match/(n_insert + n_delete + n_match + n_mismatch)
                identities.append(identity)


        step = 0.005        # bin size
        max_length = 1.0      # end of histogram

        bins = numpy.arange(0, max_length + step, step=step)

        frequencies, _ = numpy.histogram(identities, bins=bins)

        frequencies = frequencies / sum(frequencies)

        centers = (bins[:-1] + bins[1:]) / 2
        centers = centers.astype(numpy.float64)
        zeros = numpy.zeros(centers.shape)

        axes.plot(centers, frequencies, color=colors[p])
        axes.fill_between(centers, y1=zeros, y2=frequencies, color=colors[p], alpha=0.1)

        axes.set_xlabel("Identity")
        axes.set_ylabel("Proportion")

        name = path.split("/")[-1][26:].split("_VS")[0]
        print(name)
        labels.append(name)

    axes.legend(labels)

    fig.set_size_inches(12,6)
    pyplot.savefig("Read_identities.png",dpi=200)

    pyplot.show()
    pyplot.close()


def comma_separated_list(s):
    tokens = s.strip().split(",")
    return tokens

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.register("type", "comma_separated_list", comma_separated_list)

    parser.add_argument(
        "--input",
        type=comma_separated_list,
        required=True,
        help="comma separated paths of files containing read lengths"
    )
    args = parser.parse_args()

    main(paths=args.input)

"""
_VS_hg38.sorted.csv
"""