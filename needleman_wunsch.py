import argparse
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a", dest="first_sequence_file", required=True,
        help="Path to sequence file 1"
    )
    parser.add_argument(
        "-b", dest="second_sequence_file", required=True
    )

    args = parser.parse_args()

    first_sequence_file_path = os.path.abspath(args.first_sequence_file)
    second_sequence_file_path = os.path.abspath(args.second_sequence_file)

    # Extract the sequence from the sequence files
    first_sequence = extract_sequence(first_sequence_file_path)
    second_sequence = extract_sequence(second_sequence_file_path)

    print("first sequence")
    print(first_sequence)
    print(len(first_sequence))
    print("second sequence")
    print(second_sequence)
    print(len(second_sequence))


def extract_sequence(sequence_file_path):
    """
    This function takes in the path to the sequence fa file and returns the sequence as a string
    :param sequence_file_path:
    :return: sequence as a string
    """
    sequence_file = open(sequence_file_path, "r")
    # initialize the sequence string
    sequence = ""
    for line in sequence_file:
        # strip the text wrap "\n"
        line = line.rstrip()
        # skip the header line
        if line.startswith(">"):
            continue
        else:
            sequence += line
    return sequence


if __name__ == "__main__":
    main()
