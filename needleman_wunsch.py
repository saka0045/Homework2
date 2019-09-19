import argparse
import os
import numpy as np


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
    first_sequence, first_sequence_length = extract_sequence(first_sequence_file_path)
    second_sequence, second_sequence_length = extract_sequence(second_sequence_file_path)

    print("first sequence")
    print(first_sequence)
    print(first_sequence_length)
    print("second sequence")
    print(second_sequence)
    print(second_sequence_length)

    # Initialize the score matrix with zeros
    score_matrix = np.zeros([first_sequence_length + 1, second_sequence_length + 1], dtype=int)
    print(score_matrix)

    # Scoring criteria
    match = 1
    mismatch = -3
    gap_penalty = -2

    # Fill out the first column of each row with increasing gap penalty
    for row in range(0, first_sequence_length + 1):
        score_matrix[row][0] = gap_penalty * row

    # Repeat for first row in each column:
    for column in range(0, second_sequence_length + 1):
        score_matrix[0][column] = gap_penalty * column

    print(score_matrix)


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
    sequence_length = len(sequence)
    sequence_file.close()
    return sequence, sequence_length


if __name__ == "__main__":
    main()
