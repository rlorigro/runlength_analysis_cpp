from matplotlib import pyplot, colors
from datetime import datetime
import matplotlib
import argparse
import platform
import numpy
import sys
import os


if os.environ.get("DISPLAY", "") == "":
   print("no display found. Using non-interactive Agg backend")
   matplotlib.use("Agg")

if platform.system() == "Darwin":
   matplotlib.use("macosx")

A,C,G,T = 0,1,2,3

# Index key for storing base data in matrix form
BASE_TO_INDEX = {"A": 0,
                 "C": 1,
                 "G": 2,
                 "T": 3,
                 "-": 4}

INDEX_TO_BASE = ["A", "C", "G", "T"]


def get_datetime_string():
    now = datetime.now()
    now = [now.year, now.month, now.day, now.hour, now.minute, now.second, now.microsecond]
    datetime_string = "_".join(list(map(str, now)))

    return datetime_string


def normalize(frequency_matrix, pseudocount, diagonal_bias=0, logify=True):
    frequency_matrix = frequency_matrix.astype(numpy.float32)

    frequency_matrix += pseudocount

    diagonal_pseudocount = diagonal_bias
    diagonal_mask = numpy.eye(frequency_matrix.shape[0], frequency_matrix.shape[1], dtype=numpy.bool)
    frequency_matrix[diagonal_mask] += diagonal_pseudocount

    sum_y = numpy.sum(frequency_matrix, axis=1)

    probability_matrix = frequency_matrix / sum_y[:, numpy.newaxis]

    if logify:
        probability_matrix = numpy.log10(probability_matrix)

    return probability_matrix


def load_directional_base_length_matrix_from_csv(path, max_runlength_row, max_runlength_col):
    matrices = numpy.zeros([2, 4, max_runlength_row + 1, max_runlength_col + 1])
    base_index = None
    row_index = 0

    with open(path, "r") as file:
        is_data = False

        for line in file:
            if not is_data:
                if line[0] == ">":
                    # Header
                    base = line[1]
                    reversal = line.strip()[-1]

                    base_index = BASE_TO_INDEX[base]
                    reversal_index = int(reversal == "R")

                    # print(base, reversal, base_index, reversal_index)

                    is_data = True
                    row_index = 0
            else:
                if not line[0].isspace():
                    if row_index == max_runlength_row:
                        continue

                    # Data
                    row = list(map(float, line.strip().split(",")))
                    # print(base, reversal, base_index, reversal_index)
                    # print(matrices.shape)
                    matrices[reversal_index, base_index, row_index, :] = row

                    row_index += 1

                else:
                    # Space
                    is_data = False

    matrices = matrices

    return matrices


def save_directional_frequency_matrices_as_marginpolish_config(output_dir, frequency_matrices, chromosome_name=None, delimiter=",", log_normalize=False, plot=False, pseudocount=1e-12, diagonal_bias=0, default_type=int):
    if chromosome_name is not None:
        name_suffix = chromosome_name + "_"
    else:
        name_suffix = ""

    if log_normalize:
        filename = "probability_matrices_directional_MP_" + name_suffix + get_datetime_string() + ".csv"
    else:
        filename = "frequency_matrices_directional_MP_" + name_suffix + get_datetime_string() + ".csv"

    reversal_suffixes = ["F", "R"]
    output_path = os.path.join(output_dir, filename)
    file = open(output_path, "w")

    print("SAVING: %s" % output_path)

    for reversal in [0,1]:
        for base_index in range(4):
            base = INDEX_TO_BASE[base_index]
            suffix = reversal_suffixes[reversal]

            matrix = numpy.squeeze(frequency_matrices[reversal,base_index,:,:])

            type = default_type
            if log_normalize:
                matrix = normalize(matrix, pseudocount=pseudocount, diagonal_bias=diagonal_bias)
                type = float

            if plot:
                pyplot.imshow(matrix)
                pyplot.show()
                pyplot.close()

            matrix_name = "_".join([base, suffix])
            header = "\"repeatCountLogProbabilities_" + matrix_name + "\": [\n"

            file.write(header)
            row = list()
            for r in range(matrix.shape[0]):
                for x in matrix[r]:
                    row.append(str(type(x)))

            row = "    " + delimiter.join(row) + "\n],\n"

            file.write(row)

    file.close()


def main(matrix_path, output_dir, pseudocount):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print("Using pseudocount: " + str(pseudocount))

    matrix = load_directional_base_length_matrix_from_csv(path=matrix_path, max_runlength_row=50, max_runlength_col=50)

    save_directional_frequency_matrices_as_marginpolish_config(output_dir=output_dir,
                                                               frequency_matrices=matrix,
                                                               chromosome_name="genomic",
                                                               log_normalize=True,
                                                               pseudocount=pseudocount,
                                                               plot=False)


if __name__ == "__main__":
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input", "-i",
        type=str,
        required=True,
        help="Path to directional frequency matrix csv"
    )

    parser.add_argument(
        "--output", "-o",
        type=str,
        required=True,
        help="Path to output directory"
    )

    parser.add_argument(
        "--pseudocount", "-p",
        type=int,
        required=False,
        default=1,
        help="Number of pseudocounts to add to every cell in the matrix (to prevent extremely low probabilities)"
    )

    args = parser.parse_args()

    main(matrix_path=args.input, output_dir=args.output, pseudocount=args.pseudocount)