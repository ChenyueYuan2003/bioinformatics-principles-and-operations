"""
# -*- coding = utf-8 -*-
#!/usr/bin/env python
# @Project : bioinformatics
# @File : nw.py
# @Author : ycy
# @Time : 2023/5/4 9:40
# @Software : PyCharm Professional
"""

from time import time
from sys import argv


# get seq from fasta files
def read_fasta(file_path):
    sequence = ''
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            sequence += line.strip()

    return sequence


def nw(seq1, seq2):
    # Initialize scoring matrix
    n, m = len(seq1), len(seq2)
    score = [[0] * (m + 1) for _ in range(n + 1)]

    # Initialize first row and column with gap penalties
    # penalty = -8
    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] - 8
    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] - 8

    # Define BLOSUM62 matrix
    blosum62 = [[4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],
                [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
                [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],
                [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],
                [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],
                [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],
                [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],
                [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],
                [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],
                [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],
                [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],
                [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],
                [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],
                [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],
                [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2],
                [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],
                [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],
                [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],
                [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],
                [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4]]

    # The index of blosum62 matrix
    aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B',
          'Z',
          'X', '*']
    # Fill scoring matrix using BLOSUM62 scores
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = blosum62[aa.index(seq1[i - 1])][aa.index(seq2[j - 1])]
            diagonal_score = score[i - 1][j - 1] + match_score
            gap1_score = score[i - 1][j] - 8
            gap2_score = score[i][j - 1] - 8
            score[i][j] = max(diagonal_score, gap1_score, gap2_score)

    # Traceback to get alignment
    align1, align2 = '', ''
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score[i][j] == score[i - 1][j - 1] + blosum62[aa.index(seq1[i - 1])][
            aa.index(seq2[j - 1])]:
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif i > 0 and score[i][j] == score[i - 1][j] - 8:
            align1 = seq1[i - 1] + align1
            align2 = '-' + align2
            i -= 1
        else:
            align1 = '-' + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    return align1, align2, score[n][m]


def print_alignment(seq1, seq2):
    """
    Aligns two sequences by printing them out in rows of 80 characters and marks the position of each character in the entire sequence at the end of each row.
    :param seq1: The first sequence.
    :param seq2: The second sequence.
    """
    pos = 0  # The position of the current character in the entire sequence, initially 0.
    for i in range(0, len(seq1), 80):  # Processes 80 characters at a time.
        line1 = seq1[i:i + 80]  # Extracts a line of characters from the first sequence.
        line2 = seq2[i:i + 80]  # Extracts a line of characters from the second sequence.

        # Prints the first line of characters and marks the position of the current character in the entire sequence.
        for char1, char2 in zip(line1, line2):
            print(char1, end='')
            pos += 1
        print(f'  {pos:>4}')

        # Prints the second line of markers.
        for char1, char2 in zip(line1, line2):
            if char1 == char2:
                print('|', end='')
            else:
                print(' ', end='')
        print()

        # Prints the third line of characters and marks the position of the current character in the entire sequence.
        for char1, char2 in zip(line1, line2):
            print(char2, end='')
        print(f'  {pos:>4}')

        print('-' * 80)  # Prints a separator line.


def main():
    human = read_fasta(argv[1])
    mouse = read_fasta(argv[2])
    t0 = time()
    align1, align2, score = nw(human, mouse)
    print_alignment(align1, align2)
    print('score', score)
    print('time', time() - t0)


if __name__ == '__main__':
    main()
