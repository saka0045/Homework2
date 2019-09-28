"""
Script to randomly permutate the pair of sequences and align them using Needleman-Wunsch 10,000 times
Returns the alignment score for each alignment try

Discussed homework with Jeannette Rustin
"""

__author__ = "Yuta Sakai"

import argparse
import os
import random
from needleman_wunsch import extract_sequence, needleman_wunsch_algorithm


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a", dest="first_sequence_file", required=True,
        help="Path to sequence file 1"
    )
    parser.add_argument(
        "-b", dest="second_sequence_file", required=True,
        help="Path to sequence file 2"
    )
    parser.add_argument(
        "-o", dest="out_file", required=True,
        help="Path to the result output file"
    )

    args = parser.parse_args()

    # Get file path from arguments
    first_sequence_file_path = os.path.abspath(args.first_sequence_file)
    second_sequence_file_path = os.path.abspath(args.second_sequence_file)
    out_file_path = os.path.abspath(args.out_file)

    # Extract the sequence from the sequence files
    first_sequence, first_sequence_length = extract_sequence(first_sequence_file_path)
    second_sequence, second_sequence_length = extract_sequence(second_sequence_file_path)

    # Open the output result file
    out_file = open(out_file_path, "w")

    # Scoring criteria
    match = 1
    mismatch = -3
    gap_penalty = -2

    # Permutate and align the sequence 10,000 times
    n = 0

    while n in range(0, 10000):
        # Permutate the sequence
        permutated_first_sequence = permutate_sequence(first_sequence, first_sequence_length)
        permutated_second_sequence = permutate_sequence(second_sequence, second_sequence_length)

        # Get the alignment score using Needleman-Wunsch
        (alignment_score, first_aligned_sequence,
         second_aligned_sequence) = needleman_wunsch_algorithm(permutated_first_sequence, first_sequence_length,
                                                               permutated_second_sequence, second_sequence_length,
                                                               match, mismatch, gap_penalty)

        # Write the alignment score to the output file
        out_file.write(str(alignment_score) + "\n")

        n += 1

    # Close the output file
    out_file.close()
    print("Script is done running")


def permutate_sequence(sequence, sequence_length):
    """
    Function to permutate the sequence
    :param sequence:
    :param sequence_length:
    :return: randomly permutated sequence
    """
    permutated_sequence = ''.join(random.sample(sequence, sequence_length))
    return permutated_sequence


if __name__ == "__main__":
    main()
