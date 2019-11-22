from matplotlib import pyplot
from matplotlib import patches
from matplotlib.lines import Line2D
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

            else:
                coverage, error = line.strip().split(',')
                coverage = int(coverage)
                error = float(error)

                if header == "base_error":
                    base_error_x.append(coverage)
                    base_error_y.append(error)

                if header == "length_error":
                    length_error_x.append(coverage)
                    length_error_y.append(error)

    return base_error_x, base_error_y, length_error_x, length_error_y


def main(input_path):
    fig, axes = pyplot.subplots(nrows=3)

    length_match_coverage_x, length_match_coverage_y, length_mismatch_coverage_x, length_mismatch_coverage_y = read_length_coverage(
        input_path)

    base_match_coverage_x, base_match_coverage_y, base_mismatch_coverage_x, base_mismatch_coverage_y = read_base_coverage(
        input_path)

    base_error_x, base_error_y, length_error_x, length_error_y = read_percent_error(input_path)

    length_color = (0.122, 0.498, 0.584)
    base_color = (0.945, 0.71, 0.176)

    # Find sum of matches + mismatches to get total coverage
    total_coverage_y = list()
    for i in range(len(length_match_coverage_y)):
        total_coverage_y.append(length_match_coverage_y[i] + length_mismatch_coverage_y[i])

    # The lower y bound for "fill" function is just a line along y=0
    zeros = [0 for i in range(len(length_match_coverage_x))]

    # Make custom legend objects
    custom_lines = [Line2D([0], [0], color=length_color, lw=4),
                    Line2D([0], [0], color=base_color, lw=4)]

    axes[2].legend(custom_lines, ["Length", "Base"])
    axes[1].legend(custom_lines, ["Length", "Base"])

    axes[0].plot(length_match_coverage_x, total_coverage_y)
    axes[0].fill_between(length_match_coverage_x, y1=zeros, y2=length_match_coverage_y, color=length_color, alpha=0.3)

    axes[0].set_ylabel("Total Coverage")

    axes[1].plot(length_mismatch_coverage_x, length_mismatch_coverage_y)
    axes[1].fill_between(length_match_coverage_x, y1=zeros, y2=length_mismatch_coverage_y, color=length_color, alpha=0.3)

    axes[1].plot(base_mismatch_coverage_x, base_mismatch_coverage_y)
    axes[1].fill_between(base_mismatch_coverage_x, y1=zeros, y2=base_mismatch_coverage_y, color=base_color, alpha=0.3)

    axes[1].set_ylabel("Mismatches")

    axes[2].fill_between(length_error_x, y1=zeros, y2=length_error_y, color=length_color, alpha=0.3)
    axes[2].fill_between(base_error_x, y1=zeros, y2=base_error_y, color=base_color, alpha=0.3)
    axes[2].plot(length_error_x, length_error_y)
    axes[2].plot(base_error_x, base_error_y)

    axes[2].set_ylabel("P(Error)")

    axes[2].set_xlabel("Coverage")

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
