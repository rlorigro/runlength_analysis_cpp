from matplotlib import pyplot
from matplotlib import patches
import argparse


def read_length_coverage(summary_file_path):
    header = None
    length_match_coverage_x = list()
    length_match_coverage_y = list()
    length_mismatch_coverage_x = list()
    length_mismatch_coverage_y = list()

    with open(summary_file_path, 'r') as file:
        for line in file:

            if line[0] == ">":
                header = line.strip()[1:]
                print(header)

            else:
                coverage, error = line.strip().split(',')
                coverage = int(coverage)
                error = float(error)

                if header == "length_matches":
                    length_match_coverage_x.append(coverage)
                    length_match_coverage_y.append(error)

                if header == "length_mismatches":
                    length_mismatch_coverage_x.append(coverage)
                    length_mismatch_coverage_y.append(error)

    return length_match_coverage_x, length_match_coverage_y, length_mismatch_coverage_x, length_mismatch_coverage_y


def read_base_coverage(summary_file_path):
    header = None
    base_match_coverage_x = list()
    base_match_coverage_y = list()
    base_mismatch_coverage_x = list()
    base_mismatch_coverage_y = list()

    with open(summary_file_path, 'r') as file:
        for line in file:

            if line[0] == ">":
                header = line.strip()[1:]

            else:
                coverage, error = line.strip().split(',')
                coverage = int(coverage)
                error = float(error)

                if header == "base_matches":
                    base_match_coverage_x.append(coverage)
                    base_match_coverage_y.append(error)

                if header == "base_mismatches":
                    base_mismatch_coverage_x.append(coverage)
                    base_mismatch_coverage_y.append(error)

    return base_match_coverage_x, base_match_coverage_y, base_mismatch_coverage_x, base_mismatch_coverage_y


def read_percent_error(summary_file_path):
    header = None
    base_error_x = list()
    base_error_y = list()
    length_error_x = list()
    length_error_y = list()

    with open(summary_file_path, 'r') as file:
        for line in file:

            if line[0] == ">":
                header = line.strip()[1:]
                print(header)

            else:
                coverage, error = line.strip().split(',')
                coverage = int(coverage)
                error = float(error)

                print(coverage, error)

                if header == "base_error":
                    base_error_x.append(coverage)
                    base_error_y.append(error)

                if header == "length_error":
                    length_error_x.append(coverage)
                    length_error_y.append(error)

    return base_error_x, base_error_y, length_error_x, length_error_y


def main(input_path):
    fig, axes = pyplot.subplots(nrows=5)

    length_match_coverage_x, length_match_coverage_y, length_mismatch_coverage_x, length_mismatch_coverage_y = read_length_coverage(
        input_path)

    color = (0.122, 0.498, 0.584)


    zeros = [0 for i in range(len(length_match_coverage_x))]
    axes[0].plot(length_match_coverage_x, length_match_coverage_y)
    axes[0].fill_between(length_match_coverage_x, y1=zeros, y2=length_match_coverage_y, color=color, alpha=0.3)
    axes[0].set_ylabel("~Total Coverage")

    axes[1].plot(length_mismatch_coverage_x, length_mismatch_coverage_y)
    axes[1].fill_between(length_match_coverage_x, y1=zeros, y2=length_mismatch_coverage_y, color=color, alpha=0.3)
    axes[1].set_ylabel("Length Mismatch")

    base_match_coverage_x, base_match_coverage_y, base_mismatch_coverage_x, base_mismatch_coverage_y = read_base_coverage(
        input_path)

    axes[2].plot(base_mismatch_coverage_x, base_mismatch_coverage_y)
    axes[2].fill_between(base_mismatch_coverage_x, y1=zeros, y2=base_mismatch_coverage_y, color=color, alpha=0.3)
    axes[2].set_ylabel("Base Mismatch")

    base_error_x, base_error_y, length_error_x, length_error_y = read_percent_error(input_path)

    axes[3].plot(base_error_x, base_error_y)
    axes[3].fill_between(base_error_x, y1=zeros, y2=base_error_y, color=color, alpha=0.3)
    axes[3].set_ylabel("Base Error")

    axes[4].plot(length_error_x, length_error_y)
    axes[4].fill_between(length_error_x, y1=zeros, y2=length_error_y, color=color, alpha=0.3)
    axes[4].set_ylabel("Length Error")

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

    main(input_path=args.input)
