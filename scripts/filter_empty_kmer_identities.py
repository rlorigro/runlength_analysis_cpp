from matplotlib import pyplot
import argparse


def main(path):
    identities = dict()
    output_path = ".".join(path.strip().split(".")[:-1]) + "_filtered.csv"
    print("Writing file: " + output_path)

    with open(path, "r") as file, open(output_path, "w") as output_file:
        for l,line in enumerate(file):
            data = line.strip().split(",")
            kmer, n_align_match, n_insert, n_delete, n_refskip, n_softclip, n_hardclip, n_pad, n_match, n_mismatch, _ = data

            n_align_match = int(n_align_match)
            n_insert = int(n_insert)
            n_delete = int(n_delete)
            n_refskip = int(n_refskip)
            n_softclip = int(n_softclip)
            n_hardclip = int(n_hardclip)
            n_pad = int(n_pad)
            n_match = int(n_match)
            n_mismatch = int(n_mismatch)

            if n_insert + n_delete + n_match + n_mismatch == 0:
                continue
            else:
                output_file.write(line)


    return identities



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="path of file containing quadrant bounds"
    )
    args = parser.parse_args()

    main(path=args.input)