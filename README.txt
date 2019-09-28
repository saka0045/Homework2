Code to run the Needleman-Wunsch algorithm: needleman_wunsch.py
The script will take -a and -b arugments for the .fa files and if using the anchored Needleman-Wunsch, supply the
matcher text file with -t
The order of the sequence file arguments given must match the order of the columns in the matcher text file. Since the
matcher text file has Human coordinates in columns 1 and 2 and Fly coordinates in columns 3 and 4, the -a must be the
Human .fa file and -b must be the Fly .fa file, if using the matcher text:

python needleman_wunsch.py -a Human_HOX.fa -b Fly_HOX.fa -t Match_HOX.txt

Results from both HOX and PAX with non-anchored and anchored Needleman-Wunsch can be found in Aligned_Results.txt

The code to randomly permutate the pair of sequences: run_permutate.py
The pair of sequence .fa files can be supplied through the -a and -b arguments and specify the output result file with
the -o option. The result file will contain 10,000 lines with each line containing the alignment score for the randomly
permutated sequences:

python run_permutate.py -a Human_HOX.fa -b Fly_HOX.fa -o /path/to/output/file

The raw results of permutating the HOX and PAX sequences and aligning each 10,000 times are store in HOX_permutate.txt and
PAX_permutate.txt, respectively.

Histogram representing the randomly permutated alignment scores along with the original alignment score can be found in
the file HOX_Hist.pdf and PAX_Hist.pdf

For the bonus question, the aligned TITIN sequences are stored in Aligned_TITIN.txt.
