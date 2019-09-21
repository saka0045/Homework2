from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from needleman_wunsch import extract_sequence
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
    first_sequence, first_sequence_length = extract_sequence(first_sequence_file_path)
    second_sequence, second_sequence_length = extract_sequence(second_sequence_file_path)

    print("first sequence")
    print(first_sequence)
    print(first_sequence_length)
    print("second sequence")
    print(second_sequence)
    print(second_sequence_length)

    # Scoring criteria
    match = 1
    mismatch = -3
    gap_penalty = -2

    alignments = pairwise2.align.globalms(first_sequence, second_sequence, match, mismatch, gap_penalty, gap_penalty)

    for a in alignments:
        print(format_alignment(*a))


if __name__ == "__main__":
    main()
