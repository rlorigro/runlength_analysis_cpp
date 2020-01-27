import argparse
from matplotlib import pyplot

"""
EXAMPLE INPUT BLOCK:

Coverage 80
identity (M/(M+X+I+D)):	0.996161
n_matches:	5193427
n_mismatches:	3445
n_inserts:	8514
n_deletes:	8053
mismatch_rate:	0.000661
insert_rate:	0.001633
delete_rate:	0.001545
indel_rate:	0.003178
mismatches per 100kb:	66.079223
inserts per 100kb:	163.308710
deletes per 100kb:	154.466179
indels per 100kb:	317.774889

"""


def main(input_path):
    identity_per_coverage = list()

    with open(input_path, "r") as file:
        for l,line in enumerate(file):
            if line.startswith("Coverage"):
                x = int(line.strip().split(" ")[-1])
                data = list()
                data.append(x)

            elif line.startswith("mismatches per 100kb"):
                x = float(line.strip().split("\t")[-1])
                data.append(x)

            elif line.startswith("inserts per 100kb"):
                x = float(line.strip().split("\t")[-1])
                data.append(x)

            elif line.startswith("deletes per 100kb"):
                x = float(line.strip().split("\t")[-1])
                data.append(x)

            elif line == "\n":
                identity_per_coverage.append(data)

    axes = pyplot.axes()

    colors = [[0.969,0.255,0.012],
              [0.043,0.392,0.62],
              [0.427,0.800,0.108]]

    labels = ["Mismatches",
              "Inserts",
              "Deletes"]

    for i in [1,2,3]:
        x = list()
        y = list()

        for data in identity_per_coverage:
            if len(data) != 4:
                print("ERROR")
                print(data)

            x.append(data[0])
            y.append(data[i])

        axes.plot(x, y, color=colors[i-1])

    axes.legend(labels)

    axes.grid(color="k", alpha=0.1)

    axes.set_xlabel("Coverage")
    axes.set_ylabel("Error rate (per 100kb)")

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
