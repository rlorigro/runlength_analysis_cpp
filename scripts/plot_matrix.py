from matplotlib import pyplot
from matplotlib import cm
import matplotlib.colors as colors
import numpy

A,C,G,T = 0,1,2,3

# Index key for storing base data in matrix form
BASE_TO_INDEX = {"A": 0,
                 "C": 1,
                 "G": 2,
                 "T": 3,
                 "-": 4}

INDEX_TO_BASE = ["A", "C", "G", "T"]


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


def truncate_colormap(cmap, min=0.0, max=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=min, b=max),
        cmap(numpy.linspace(min, max, n)))
    return new_cmap


def calculate_diagonal_accuracy(matrix, axis=1):
    sums = numpy.sum(matrix, axis=axis)
    sums = numpy.expand_dims(sums, axis)

    diagonal_mask = numpy.eye(matrix.shape[0], matrix.shape[1], dtype=numpy.bool)

    diagonals = matrix[diagonal_mask]
    diagonals = numpy.expand_dims(diagonals, axis)

    accuracies = diagonals/sums

    # print(sums)
    # print(sums.shape)
    # print("accuracies on axis", axis)
    # print(accuracies)
    # print(accuracies.shape)

    total_accuracy = numpy.sum(diagonals)/numpy.sum(matrix)
    print(total_accuracy)

    return accuracies


def calculate_average_diagonal_distance(matrix, axis=1):
    sums = numpy.sum(matrix, axis=axis)
    sums = numpy.expand_dims(sums, axis)

    weight_mask = numpy.zeros(matrix.shape)
    for i0 in range(weight_mask.shape[0]):
        for i1 in range(weight_mask.shape[1]):
            weight_mask[i0][i1] = abs(i1-i0)

    if axis == 1:
        matrix = matrix/sums
    elif axis == 0:
        matrix = matrix/sums

    # test_figure = pyplot.figure()
    # test_axes = pyplot.axes()
    # test_axes.imshow(matrix)
    # pyplot.show()
    # pyplot.close()

    total_off_diagonal_distance = numpy.sum(weight_mask*matrix)

    return


def main():

    # paths = [
    #     "/home/ryan/data/lc2020/length_frequency_matrix_nondirectional_shasta_r9_305.csv",
    #     "/home/ryan/data/lc2020/length_frequency_matrix_nondirectional_shasta_r9_360.csv",
    #     "/home/ryan/data/lc2020/length_frequency_matrix_nondirectional_shasta_r10_348.csv"
    # ]
    #
    # labels = [
    #     "Guppy v3.0.5 r9.4",
    #     "Guppy v3.6.0 r9.4",
    #     "Guppy v3.4.8 r10.3"
    # ]

    # ------------------------------------------------------------------------------------------------------

    # paths = ["/home/ryan/data/lc2020/length_frequency_matrix_nondirectional_shasta_r9_360.csv",
    #          "/home/ryan/data/lc2020/length_frequency_matrix_nondirectional_shasta_r10_348.csv"]
    #
    # labels = ["R9.4.1 Guppy v3.6.0",
    #           "R10.3 Guppy v3.4.8"]

    # ------------------------------------------------------------------------------------------------------

    # paths = [
    #     # "/home/ryan/data/lc2020/asm_rle/shasta_HG002_assembly_length_frequency_matrices_guppy_235_r9.csv",
    #     "/home/ryan/data/lc2020/asm_rle/shasta_HG002_assembly_length_frequency_matrices_guppy_305_r9.csv",
    #     "/home/ryan/data/lc2020/asm_rle/shasta_HG002_assembly_length_frequency_matrices_guppy_360_r9.csv",
    #     # "/home/ryan/data/lc2020/asm_rle/shasta_HG002_assembly_length_frequency_matrices_guppy_348_r10.csv"
    #     ]
    #
    # labels = [
    #     # "Guppy v2.3.5 r9.4",
    #     "Guppy v3.0.5 r9.4",
    #     "Guppy v3.6.0 r9.4",
    #     # "Guppy v3.4.8 r10.3"
    # ]

    # ------------------------------------------------------------------------------------------------------

    paths = [
        "/home/ryan/code/runlength_analysis_cpp/build/output/rle_test/length_frequency_matrix_nondirectional.csv",
        "/home/ryan/code/runlength_analysis_cpp/build/output/rle_test/length_frequency_matrix_nondirectional.csv",
        ]

    labels = [
        "test",
        "test",
    ]


    figure,axes = pyplot.subplots(ncols=len(paths), nrows=1)
    colormap = pyplot.get_cmap("viridis")
    # colormap = truncate_colormap(colormap, max=0.5)

    # colormap = pyplot.get_cmap("twilight")
    # colormap = truncate_colormap(colormap, max=0.5)

    accuracies_per_sample = list()

    for p,path in enumerate(paths):
        matrix = load_base_length_matrix_from_csv(path,
                                                  max_runlength_row=50,
                                                  max_runlength_col=50)

        matrix = numpy.sum(matrix, axis=0)[1:31,1:31]
        accuracies = calculate_diagonal_accuracy(matrix,axis=1)
        accuracies_per_sample.append(accuracies)
        calculate_average_diagonal_distance(matrix)

        sums = numpy.sum(matrix, axis=1)
        sums = numpy.expand_dims(sums, 1)

        matrix = matrix/sums
        # matrix += 1

        axes[p].imshow(matrix, cmap=colormap)

        axes[p].set_xticks([0,5-1,10-1,15-1,20-1,25-1,30-1]) #,30-1,35-1])
        axes[p].set_xticklabels([1,5,10,15,20,25,30]) #,30,35])
        axes[p].set_yticks([0,5-1,10-1,15-1,20-1,25-1,30-1]) #,30-1,35-1])
        axes[p].set_yticklabels([1,5,10,15,20,25,30]) #,30,35])

        axes[p].set_xlabel("Predicted Length (bp)", fontsize=16)
        axes[p].set_title(labels[p], fontsize=22)

    axes[0].set_ylabel("True Length (bp)", fontsize=16)
    figure.set_size_inches(len(paths)*8, 8)
    pyplot.savefig("rle.pdf")
    pyplot.savefig("rle.png", dpi=200)

    pyplot.show()
    pyplot.close()

    figure2 = pyplot.figure()
    axes2 = pyplot.axes()

    for accuracies in accuracies_per_sample:
        axes2.plot(accuracies)

    axes2.legend(labels)
    axes2.set_xlabel("Homopolymer length")
    axes2.set_ylabel("Proportion Correct")
    figure2.set_size_inches(8,6)
    pyplot.savefig("rle_accuracy.png", dpi=200)
    pyplot.show()
    pyplot.close()


if __name__ == "__main__":
    main()
