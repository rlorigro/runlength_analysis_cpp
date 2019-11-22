from os.path import isfile, isdir, join, dirname, exists, basename
from os import listdir, remove, walk, mkdir
from matplotlib import pyplot
from matplotlib import patches
import argparse
import numpy
import math
from matplotlib import pyplot
import numpy


def weibull_cdf(x, shape, scale):
    """  Tail cumulative probability of Weibull distribution

        :param x: points at which to evaluate
        :param shape: Shape parameter
        :param scale: Scale parameter
    """
    return numpy.exp(-numpy.power(x / scale, shape))


def evaluate_discrete_weibull(scale, shape, x):
    """
    Evaluate the discrete form of the distribution at points specified by x vector
    :param scale:
    :param shape:
    :param x:
    :return:
    """
    distribution = weibull_cdf(x, shape, scale) - weibull_cdf(x + 1, shape, scale)

    return distribution


def calculate_mode(scale, shape):
    """
    From wiki, lambda=scale, k=shape
    :param scale:
    :param shape:
    :param x:
    :return:
    """

    if shape > 1:
        mode = scale * ((shape-1) / shape)**(1 / shape)
        mode = max(0, mode)

    else:
        mode = 0

    return mode


def main():
    pass


if __name__ == "__main__":
    main()


def calculate_entropy(pdf):
    pdf += 1e-6

    log_pdf = numpy.log(pdf)

    h_max = math.log(len(pdf))

    h = numpy.sum((log_pdf*pdf)) / h_max

    return h


def get_all_file_paths_by_type(parent_directory_path, file_extension, sort=True, recursive=True):
    """
    Given a parent directory, iterate all files within, and return those that end in the extension provided by user.
    File paths returned in sorted order by default.
    :param parent_directory_path:
    :param file_extension:
    :param sort:
    :return:
    """
    all_files = set()

    for root, dirs, files in walk(parent_directory_path):
        for file in files:
            if file.endswith(file_extension):
                all_files.add(join(root, file))

        if not recursive:
            break

    all_files = list(all_files)

    if sort:
        all_files.sort()

    return all_files


def read_bounds(input_path):
    centers = list()
    half_widths = list()

    with open(input_path, 'r') as file:
        for line in file:
            center, half_size = line.strip().split('\t')

            half_size = float(half_size)

            center = center.split(',')
            center = list(map(float, center))

            centers.append(center)
            half_widths.append(half_size)

    return centers, half_widths


def main(input_path, axes):
    centers, half_widths = read_bounds(input_path)

    # axes = pyplot.axes()

    color_map = pyplot.get_cmap('viridis')

    y_max = 0
    x_max = 0
    x = numpy.arange(0, 50)
    for i in range(len(centers)):
        corner = (centers[i][0] - half_widths[i], centers[i][1] - half_widths[i])
        scale = centers[i][0]
        shape = centers[i][1]
        width = half_widths[i]*2

        y = evaluate_discrete_weibull(shape=shape, scale=scale, x=x)
        entropy = calculate_entropy(y)
        # print(entropy)
        color = color_map(-entropy)

        rectangle = patches.Rectangle(corner, width, width, linewidth=0, edgecolor='r', facecolor=color)
        axes.add_patch(rectangle)

        right = centers[i][0] + half_widths[i]
        top = centers[i][1] + half_widths[i]

        if right > x_max:
            x_max = right

        if top > y_max:
            y_max = top

    axes.set_xlim(0, x_max)
    axes.set_ylim(0, y_max)

    # pyplot.show()
    # pyplot.close()


def plot_all(input_dir):
    paths = get_all_file_paths_by_type(parent_directory_path=input_dir, file_extension=".txt", sort=False, recursive=False)

    paths = [path for path in paths if basename(path).startswith("quad_tree_")]

    paths = sorted(paths, key=lambda x: int(x.split(".")[0].split("_")[-1]))

    # print(paths)

    fig = pyplot.figure()
    axes = pyplot.axes()
    pyplot.ion()

    for p,path in enumerate(paths):
        main(input_path=path, axes=axes)

        axes.set_xlim([0, 14])
        axes.set_ylim([0, 30])
        pyplot.draw()
        pyplot.pause(0.0001)

        if (p == len(paths)-1):
            fig.savefig(join(input_dir, "quad_compression.png"), dpi=200)
            # pyplot.show()

        pyplot.cla()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path of DIRECTORY containing quadrant bounds with prefix 'quad_tree_` and suffix '.txt'"
    )
    args = parser.parse_args()

    plot_all(input_dir=args.input)
