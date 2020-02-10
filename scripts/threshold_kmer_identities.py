from matplotlib import pyplot
import argparse


def main(path, identity_cutoff, match_cutoff):
    output_stats_path = ".".join(path.strip().split(".")[:-1]) + "_identity-" + str(identity_cutoff) + "_match-" + str(match_cutoff) + ".csv"
    output_kmers_path = ".".join(path.strip().split(".")[:-1]) + "_identity-" + str(identity_cutoff) + "_match-" + str(match_cutoff) + "_KMERS_ONLY.csv"

    print("Writing file: " + output_stats_path)

    n_found = 0
    with open(path, "r") as file, open(output_stats_path, "w") as output_stats_file, open(output_kmers_path, "w") as output_kmers_file:
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

            identity = float(n_match)/(n_insert + n_delete + n_match + n_mismatch)

            if identity > identity_cutoff and n_match < match_cutoff:
                output_stats_file.write(line.strip())
                output_stats_file.write(str(identity))
                output_stats_file.write("\n")

                output_kmers_file.write(kmer)
                output_kmers_file.write("\n")

                n_found += 1

        print(n_found)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="path of file containing quadrant bounds"
    )
    parser.add_argument(
        "--identity",
        type=float,
        required=True,
        help="Decimal cutoff between 0 and 1 for identity"
    )
    parser.add_argument(
        "--matches",
        type=int,
        required=False,
        default=0,
        help="Integer cutoff for total matches found in the kmer"
    )
    args = parser.parse_args()

    main(path=args.input, identity_cutoff=args.identity, match_cutoff=args.matches)