from matplotlib import pyplot
from matplotlib import patches
import argparse


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


def main(input_path):
    centers, half_widths = read_bounds(input_path)

    axes = pyplot.axes()

    y_max = 0
    x_max = 0
    for i in range(len(centers)):
        corner = (centers[i][0] - half_widths[i], centers[i][1] - half_widths[i])
        width = half_widths[i]*2

        rectangle = patches.Rectangle(corner, width, width, linewidth=1, edgecolor='r', facecolor='none')
        axes.add_patch(rectangle)

        right = centers[i][0] + half_widths[i]
        top = centers[i][1] + half_widths[i]

        if right > x_max:
            x_max = right

        if top > y_max:
            y_max = top

    axes.set_xlim(0, x_max)
    axes.set_ylim(0, y_max)

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
