from matplotlib import pyplot, colors
from datetime import datetime
import matplotlib
import argparse
import platform
import numpy
import sys
import os


HG38_PRIOR_STRING = ">AT prior\n" + \
                    "-9.064598218,-0.156115265,-0.733150752,-1.137260069,-1.562448179,-1.972014127,-2.496940054,-2.876194859,-3.319246222,-3.562261244,-3.751292148,-3.966307608,-4.075210701,-4.129711096,-4.180840284,-4.230209309,-4.296352417,-4.390849486,-4.486302913,-4.558093185,-4.630189010,-4.682977808,-4.720874604,-4.777244445,-4.840324204,-4.925033952,-5.033877089,-5.142547816,-5.237552201,-5.351023680,-5.484129434,-5.681321567,-5.827054480,-5.950654866,-6.085505317,-6.128587422,-6.130605054,-6.121103702,-6.162595327,-6.235938321,-6.322659140,-6.480266994,-6.630029314,-6.761402160,-6.857772342,-6.903230216,-6.897280883,-6.954008508,-7.031174462,-7.183784626,-7.265257668\n\n" + \
                    ">GC prior\n" + \
                    "-8.952736661,-0.127049254,-0.724514025,-1.304860743,-1.921083364,-2.548855332,-3.251059900,-4.067296280,-4.767137850,-5.178365701,-5.397884227,-5.602488643,-5.813172395,-6.043180632,-6.276043051,-6.532780912,-6.740549057,-6.903518638,-7.126661858,-7.271495424,-7.554796652,-7.630517366,-8.049646674,-7.998494152,-8.049646674,-8.475615406,-8.952736661,-8.174585411,-8.952736661,-8.174585411,-8.651706665,-8.952736661,-8.952736661,-8.952736661,-8.952736661,-8.952736661,-8.952736661,-8.952736661,-8.952736661,-8.952736661,-8.952736661,-8.952736661,-8.952736661,-8.952736661,-8.952736661,-8.952736661,-8.952736661,-8.952736661,-8.952736661,-8.952736661,-8.952736661\n\n"

ECOLI_K12_PRIOR_STRING = ">AT prior\n" + \
                    "-6.205901275,-0.146046406,-0.713553033,-1.204055110,-1.715900633,-2.136270172,-2.623837912,-3.229551296,-3.892034055,-4.950628770,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275,-6.205901275\n\n" + \
                    ">GC prior\n" + \
                    "-6.258995296,-0.122685760,-0.691951276,-1.458793905,-2.176425332,-2.928175830,-3.771856921,-4.510807269,-5.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296,-6.258995296\n\n"


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


def load_base_length_matrix_from_csv(path, max_runlength_row, max_runlength_col):
    matrices = numpy.zeros([4, max_runlength_row + 1, max_runlength_col + 1])
    base_index = None
    row_index = 0

    with open(path, "r") as file:
        is_data = False

        for line in file:
            if not is_data:
                if line[0] == ">":
                    # Header
                    base = line[1]
                    base_index = BASE_TO_INDEX[base]

                    is_data = True
                    row_index = 0
            else:
                if not line[0].isspace():
                    if row_index == max_runlength_row:
                        continue

                    # Data
                    row = list(map(float, line.strip().split(",")))

                    matrices[base_index, row_index, :] = row

                    row_index += 1


                else:
                    # Space
                    is_data = False

    matrices = matrices

    return matrices


def save_nondirectional_frequency_matrices_as_delimited_text(output_dir, prior, name, frequency_matrices, chromosome_name=None, delimiter=",", log_normalize=False, pseudocount=1e-12, diagonal_bias=0, plot=False, default_type=int, filename=None):
    if filename is None:
        if chromosome_name is not None:
            name_suffix = chromosome_name + "_"
        else:
            name_suffix = ""

        if log_normalize:
            filename = "probability_matrices_" + name_suffix + get_datetime_string() + ".csv"
        else:
            filename = "frequency_matrices_" + name_suffix + get_datetime_string() + ".csv"

    output_path = os.path.join(output_dir, filename)
    file = open(output_path, "w")

    file.write(">Name\n")
    file.write(name + " with pseudocounts " + str(pseudocount) + "\n\n")

    if prior == "ecoli":
        file.write(ECOLI_K12_PRIOR_STRING)
    elif prior == "human":
        file.write(HG38_PRIOR_STRING)
    else:
        exit("ERROR: unknown or unspecified prior: " + prior)

    for base_index in range(4):
        base = INDEX_TO_BASE[base_index]

        matrix = numpy.squeeze(frequency_matrices[base_index,:,:])

        type = default_type
        if log_normalize:
            matrix = normalize(matrix, pseudocount=pseudocount, diagonal_bias=diagonal_bias)
            type = float

        if plot:
            pyplot.imshow(matrix)
            pyplot.show()
            pyplot.close()

        matrix_name = base

        if log_normalize:
            matrix_name += " likelihood"

        header = ">" + matrix_name + "\n"

        file.write(header)

        for r in range(matrix.shape[0]):
            row = [str(type(x)) for x in matrix[r]]
            row = delimiter.join(row) + "\n"

            file.write(row)

        file.write("\n")

    file.close()


def main(matrix_path, output_dir, pseudocount, prior, name):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print("Using pseudocount: " + str(pseudocount))

    matrix = load_base_length_matrix_from_csv(path=matrix_path, max_runlength_row=50, max_runlength_col=50)

    output_filename_prefix = matrix_path.split("/")[-1]
    output_filename_prefix = ".".join(output_filename_prefix.split(".")[:-1])
    filename = output_filename_prefix + "_shasta_bayesian_config.csv"

    print(filename)

    save_nondirectional_frequency_matrices_as_delimited_text(output_dir=output_dir,
                                                             frequency_matrices=matrix,
                                                             chromosome_name="genomic",
                                                             log_normalize=True,
                                                             filename=filename,
                                                             pseudocount=pseudocount,
                                                             plot=False,
                                                             prior=prior,
                                                             name=name)




if __name__ == "__main__":
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--prior",
        type=str,
        required=True,
        help="Choose between Ecoli ('ecoli') and Human ('human') priors for runlength"
    )
    parser.add_argument(
        "--name", "-n",
        type=str,
        required=True,
        help="A descriptive name to be put in the header of the config file. Ideally describing the organism + region + sequencer"
    )
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

    main(matrix_path=args.input, output_dir=args.output, pseudocount=args.pseudocount, prior=args.prior, name=args.name)