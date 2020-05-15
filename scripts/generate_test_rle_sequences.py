from random import randint


def complement(base):
    if base == "A":
        return "T"
    elif base == "C":
        return "G"
    elif base == "G":
        return "C"
    elif base == "T":
        return "A"
    else:
        exit("ERROR: unidentified base has no complement: " + base)


def choose_next_base(base):
    if base == "A":
        return a_transitions[randint(0,2)]
    elif base == "C":
        return c_transitions[randint(0,2)]
    elif base == "G":
        return g_transitions[randint(0,2)]
    elif base == "T":
        return t_transitions[randint(0,2)]
    else:
        exit("ERROR: unidentified base has no complement: " + base)


def choose_sequence(length):
    # Start with random base
    sequence = bases[randint(0,3)]

    # Given the current base choose a random base from the set of other bases
    for i in range(length-1):
        sequence += choose_next_base(sequence[-1])

    return sequence


def write_reference_sequence_to_file(reference_sequence):
    with open("test_rle_reference.fasta", "w") as file:
        file.write(">0\n")
        file.write("".join(reference_sequence))
        file.write("\n")


def write_sequences_to_file(sequences):
    with open("test_rle_sequences.fasta", "w") as file:
        for i in range(len(sequences)):
            file.write(">")
            file.write(str(i+1))
            file.write("\n")
            file.write("".join(sequences[i]))
            file.write("\n")


def main():
    n_samples = 20
    max_length = 50
    n_bases = max_length*max_length*n_samples

    reference_bases = choose_sequence(n_bases)
    sequences = [list() for i in range(max_length)]
    reference_sequence = list()

    for i_base in range(n_bases):
        reference_sequence.append(reference_bases[i_base]*((i_base % max_length)+1))

    for i_seq in range(max_length):
        for i_base in range(n_bases):
            sequences[i_seq].append(reference_bases[i_base]*(i_seq+1))

    write_sequences_to_file(sequences)
    write_reference_sequence_to_file(reference_sequence)

if __name__ == "__main__":
    bases = ["A","C","G","T"]

    a_transitions = ["C","G","T"]
    c_transitions = ["A","G","T"]
    g_transitions = ["C","A","T"]
    t_transitions = ["C","G","A"]

    main()
