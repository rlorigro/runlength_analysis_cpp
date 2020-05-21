from matplotlib import pyplot
import argparse
import numpy
import os


def find_coverage_intervals(lengths, genome_size=3200000000):
    length_sum = 0
    total_coverage = sum(lengths)/genome_size
    print(total_coverage)

    coverages = list(range(int(round(total_coverage/10))+2))
    coverages = [c*10 for c in coverages]

    print(coverages)

    for length in reversed(lengths):
        length_sum += length
        if (length_sum/genome_size) > coverages[0]:
            print(str(coverages[0]) + "x coverage greater than " + str(length))

            del coverages[0]

            if len(coverages) == 0:
                break


def plot_histogram(axes, lengths, i):
    colors = [(0.945, 0.267, 0.176), (0.122, 0.498, 0.584), (0.945, 0.71, 0.176)]

    step = 1000        # bin size
    max_length = 300000      # end of histogram

    bins = numpy.arange(0, max_length + step, step=step)

    frequencies, _ = numpy.histogram(lengths, bins=bins)


    centers = (bins[:-1] + bins[1:]) / 2
    centers = centers.astype(numpy.int64)
    zeros = numpy.zeros(centers.shape)

    frequencies *= centers

    axes.plot(centers, frequencies, color=colors[i])
    axes.fill_between(centers, y1=zeros, y2=frequencies, color=colors[i], alpha=0.3)

    axes.set_xlabel("Length (bp)")
    axes.set_ylabel("Coverage (bp)")

    # pyplot.show()
    # pyplot.close()


def main(paths):
    axes = pyplot.axes()

    for p,path in enumerate(paths):
        print(path)

        with open(path, "r") as file:
            for l,line in enumerate(file):
                if l == 1:
                    lengths = list(map(int,line.strip().split(",")))
                    lengths = [l for l in lengths if l > 100]
                    find_coverage_intervals(lengths)

            plot_histogram(axes=axes, lengths=lengths, i=p)



    axes.legend(["PAD65442", "PAD64459", "PAD64591"])
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

