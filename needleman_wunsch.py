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
        "-b", dest="second_sequence_file", required=True,
        help="Path to sequence file 2"
    )
    parser.add_argument(
        "-t", dest="text_file", required=False,
        help="Path to the matchers.txt, having this argument will enable the anchored global alignment"
    )

    args = parser.parse_args()

    first_sequence_file_path = os.path.abspath(args.first_sequence_file)
    second_sequence_file_path = os.path.abspath(args.second_sequence_file)

    # Extract the sequence from the sequence files
    first_sequence, first_sequence_length = extract_sequence(first_sequence_file_path)
    second_sequence, second_sequence_length = extract_sequence(second_sequence_file_path)

    # Scoring criteria
    match = 1
    mismatch = -3
    gap_penalty = -2

    if args.text_file is not None:
        print("Performing anchored global alignment")
        text_file_path = os.path.abspath(args.text_file)
        print("Using text file to retrieve start and stop position at: " + text_file_path)

        print("Aligning the following sequences using anchored Needleman-Wunsch algorithm")
        print("First sequence:")
        print(first_sequence)
        print("Second sequence:")
        print(second_sequence)

        text_file = open(text_file_path, "r")
        anchored_alignement_score, anchored_first_aligned_sequence, anchored_second_aligned_sequence = anchored_needleman_wunsch(
            first_sequence, first_sequence_length, second_sequence, second_sequence_length, text_file, match, mismatch,
        gap_penalty)

        print("Alignment score is :" + str(anchored_alignement_score))
        print("Aligned sequences using anchored Needleman-Wunsch:\n")
        print(anchored_first_aligned_sequence + "\n")
        print(anchored_second_aligned_sequence)

        text_file.close()

    # If no -t argument, proceed with non-anchored Needleman-Wunsch
    else:
        print("Aligning the following sequences using Needleman-Wunsch algorithm")
        print("First sequence:")
        print(first_sequence)
        print("Second sequence:")
        print(second_sequence)

        # Align using non-anchored version of Needleman-Wunsch
        alignment_score, first_aligned_sequence, second_aligned_sequence = needleman_wunsch_algorithm(first_sequence,
                                                                                                      first_sequence_length,
                                                                                                      second_sequence,
                                                                                                      second_sequence_length,
                                                                                                      match, mismatch,
                                                                                                      gap_penalty)

        print("Alignment score is :" + str(alignment_score))
        print("Aligned sequences using non-anchored Needleman-Wunsch:\n")
        print(first_aligned_sequence + "\n")
        print(second_aligned_sequence)


def anchored_needleman_wunsch(first_sequence, first_sequence_length, second_sequence, second_sequence_length,
                              text_file, match, mismatch, gap_penalty):
    """
    Performs the anchored Needleman-Wunsch, requires -t argument when invoking this script
    :param first_sequence:
    :param first_sequence_length:
    :param second_sequence:
    :param second_sequence_length:
    :param text_file:
    :return: Returns the alignment score and the sequences using anchored Needleman-Wunsch
    """
    # Store the start and stop positions in a list
    list_of_coordinates = []
    for line in text_file:
        # Remove any "\n" at the end of line
        line = line.rstrip()
        coordinates_per_line = line.split("\t")
        # Convert the coordinates to integers
        coordinates_per_line = list(map(int, coordinates_per_line))
        list_of_coordinates.append(coordinates_per_line)
    print("Coordinates of matched regions:")
    print(list_of_coordinates)
    # Initialize for anchored global alignment
    first_sequence_starting_position = 0
    second_sequence_starting_position = 0
    anchored_alignement_score = 0
    anchored_first_aligned_sequence = ""
    anchored_second_aligned_sequence = ""
    for coordinates in list_of_coordinates:
        # Align the sequences before the set of coordinates with Needleman-Wunsch
        first_sequence_to_align = first_sequence[first_sequence_starting_position:coordinates[0] - 1]
        second_sequence_to_align = second_sequence[second_sequence_starting_position:coordinates[2] - 1]
        print("Aligning the following sequences using Needleman-Wunsch:")
        print(first_sequence_to_align)
        print(second_sequence_to_align + "\n")
        (alignment_score, first_aligned_sequence,
         second_aligned_sequence) = needleman_wunsch_algorithm(first_sequence_to_align,len(first_sequence_to_align),
                                                               second_sequence_to_align,len(second_sequence_to_align),
                                                               match, mismatch, gap_penalty)
        print("Aligned sequences:")
        print(first_aligned_sequence)
        print(second_aligned_sequence + "\n")

        # Add the alignment score and extend the aligned sequence
        anchored_alignement_score += alignment_score
        anchored_first_aligned_sequence += first_aligned_sequence
        anchored_second_aligned_sequence += second_aligned_sequence

        # Extract the matched region from the sequence
        first_matched_region = first_sequence[coordinates[0] - 1:coordinates[1]]
        second_matched_region = second_sequence[coordinates[2] - 1:coordinates[3]]
        print("Extracting the following matched regions")
        print(first_matched_region)
        print(second_matched_region + "\n")
        # Calculate the alignment score for the matched region
        anchored_alignement_score += len(first_matched_region) * match
        # Add the matched region to the aligned sequence
        anchored_first_aligned_sequence += first_matched_region
        anchored_second_aligned_sequence += second_matched_region

        # Move the starting positions to so the next loop can start where the matched region left off
        first_sequence_starting_position = coordinates[1]
        second_sequence_starting_position = coordinates[3]
    # Align the rest of the sequence using Needleman-Wunsch
    rest_of_first_sequence = first_sequence[first_sequence_starting_position:first_sequence_length]
    rest_of_second_sequence = second_sequence[second_sequence_starting_position:second_sequence_length]
    print("Aligning the rest of the sequence using Needleman-Wunsch:")
    print(rest_of_first_sequence)
    print(rest_of_second_sequence + "\n")
    (last_alignment_score,
     last_first_alignment_sequence,
     last_second_alignment_sequence) = needleman_wunsch_algorithm(rest_of_first_sequence,len(rest_of_first_sequence),
                                                                  rest_of_second_sequence, len(rest_of_second_sequence),
                                                                  match, mismatch, gap_penalty)
    print("Aligned sequence:")
    print(last_first_alignment_sequence)
    print(last_second_alignment_sequence + "\n")
    # Calculate the alignment score of the last region and add the sequences to the aligned sequence
    anchored_alignement_score += last_alignment_score
    anchored_first_aligned_sequence += last_first_alignment_sequence
    anchored_second_aligned_sequence += last_second_alignment_sequence
    return anchored_alignement_score, anchored_first_aligned_sequence, anchored_second_aligned_sequence


def needleman_wunsch_algorithm(first_sequence, first_sequence_length, second_sequence, second_sequence_length,
                               match, mismatch, gap_penalty):
    """
    Performs global alignment using the Needleman-Wunsch algorithm
    :param first_sequence:
    :param first_sequence_length:
    :param second_sequence:
    :param second_sequence_length:
    :return: Aligned sequences and the alignment score
    """
    # Initialize the score matrix with zeros
    score_matrix = np.zeros([first_sequence_length + 1, second_sequence_length + 1], dtype=int)
    traceback_matrix = np.empty([first_sequence_length + 1, second_sequence_length + 1], dtype=object)

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
    # Get the alignment score
    alignment_score = score_matrix[first_sequence_length][second_sequence_length]
    # print(score_matrix)
    # print(traceback_matrix)
    # Traceback
    # Set and initialize the traceback conditions
    i = first_sequence_length
    j = second_sequence_length
    first_aligned_sequence = ""
    second_aligned_sequence = ""
    # Start at the bottom right corner of the traceback_matrix
    # While loop before it reaches any of the sides or the top left corner
    while i > 0 and j > 0:
        # If the direction is diagonal, add the corresponding sequence from both sequences
        if traceback_matrix[i][j] == "diagonal":
            first_aligned_sequence += first_sequence[i - 1]
            second_aligned_sequence += second_sequence[j - 1]
            i -= 1
            j -= 1
        # If the direction is left, add a gap to the first sequence (on the side of table)
        # but take the sequence from the second sequence (on the top of the table)
        elif traceback_matrix[i][j] == "left":
            first_aligned_sequence += "-"
            second_aligned_sequence += second_sequence[j - 1]
            j -= 1
        # If the direction is up, take the sequence from the first sequence and add a gap to the second sequence
        elif traceback_matrix[i][j] == "up":
            first_aligned_sequence += first_sequence[i - 1]
            second_aligned_sequence += "-"
            i -= 1
    # If it reached the left end of the table, add the corresponding sequence to the first sequence and add a gap
    # to the second sequence and move up
    while i > 0:
        first_aligned_sequence += first_sequence[i - 1]
        second_aligned_sequence += "-"
        i -= 1
    # If it reached the top of the table, add a gap to the first sequence and add the corresponding sequence from the
    # second sequence and move left
    while j > 0:
        first_aligned_sequence += "-"
        second_aligned_sequence += second_sequence[j - 1]
        j -= 1
    # Invert the sequences since we started at the end
    first_aligned_sequence = first_aligned_sequence[::-1]
    second_aligned_sequence = second_aligned_sequence[::-1]
    return alignment_score, first_aligned_sequence, second_aligned_sequence


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
    #print("out of " + str(diagonal_score) + ", " + str(up_score) + ", " + str(left_score) + ", max is " + str(score))
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
