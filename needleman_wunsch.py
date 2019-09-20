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
    traceback_matrix = np.empty([first_sequence_length + 1, second_sequence_length + 1], dtype=object)

    # Scoring criteria
    match = 1
    mismatch = -3
    gap_penalty = -2

    # Fill out the first column of each row with increasing gap penalty for score matrix
    # Fill out the first column of each row with "up" for traceback matrix
    for row in range(0, first_sequence_length + 1):
        score_matrix[row][0] = gap_penalty * row
        traceback_matrix[row][0] = "up"

    # Repeat for first row in each column for score matrix and fill out "left" for traceback matrix
    for column in range(0, second_sequence_length + 1):
        score_matrix[0][column] = gap_penalty * column
        traceback_matrix[0][column] = "left"

    # Fill out the first cell for traceback matrix
    traceback_matrix[0][0] = "start"

    # Calculate the score for the gut of the matrix
    for row in range(1, first_sequence_length + 1):
        for column in range(1, second_sequence_length + 1):
            score, direction = calculate_matrix_score(first_sequence, row, second_sequence, column, match,
                                                      mismatch, gap_penalty, score_matrix)
            score_matrix[row][column] = score
            traceback_matrix[row][column] = direction

    print(score_matrix)
    print(traceback_matrix)

    # Traceback
    i = first_sequence_length
    j = second_sequence_length
    first_aligned_sequence = ""
    second_aligned_sequence = ""

    while i > 0 and j > 0:
        if traceback_matrix[i][j] == "diagonal":
            first_aligned_sequence += first_sequence[i - 1]
            second_aligned_sequence += second_sequence[j - 1]
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == "left":
            first_aligned_sequence += "-"
            second_aligned_sequence += second_sequence[j - 1]
            j -= 1
        elif traceback_matrix[i][j] == "up":
            first_aligned_sequence += first_sequence[i - 1]
            second_aligned_sequence += "-"
            i -= 1

    while i > 0:
        first_aligned_sequence += first_sequence[i - 1]
        second_aligned_sequence += "-"
        i -= 1

    while j > 0:
        first_aligned_sequence += "-"
        second_aligned_sequence += second_sequence[j - 1]
        j -= 1

    first_aligned_sequence = first_aligned_sequence[::-1]
    second_aligned_sequence = second_aligned_sequence[::-1]

    print(first_aligned_sequence)
    print(second_aligned_sequence)


def calculate_matrix_score(first_sequence, row, second_sequence, column, match, mismatch, gap_penalty, score_matrix):
    """
    Function to calculate the matrix used for Needleman-Wunsch algorithm. Returns the highest value from either
    coming to the cell diagonally, from the top or from the left cell
    :param first_sequence:
    :param row:
    :param second_sequence:
    :param column:
    :param match:
    :param mismatch:
    :param gap_penalty:
    :param score_matrix:
    :return: highest score for the given cell in the matrix
    """
    # To calculate the diagonal score, see if the sequence match up
    # If the sequences match, add the match score to the adjacent top left cell
    if first_sequence[row - 1] == second_sequence[column - 1]:
        diagonal_score = score_matrix[row - 1][column - 1] + match
    # If the sequences don't match, add the mismatch score
    else:
        diagonal_score = score_matrix[row - 1][column - 1] + mismatch
    # To calculate the score coming from the left, add the gap penalty to the left cell
    left_score = score_matrix[row][column - 1] + gap_penalty
    # To calculate the score coming from the top, add the gap penalty to the top cell
    up_score = score_matrix[row - 1][column] + gap_penalty
    # Assess which score is the max
    score = max(diagonal_score, left_score, up_score)
    if diagonal_score == score:
        direction = "diagonal"
    elif left_score == score:
        direction = "left"
    elif up_score == score:
        direction = "up"
    print("out of " + str(diagonal_score) + ", " + str(up_score) + ", " + str(left_score) + ", max is " + str(score))
    return score, direction


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
