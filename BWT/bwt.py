"""
# -*- coding = utf-8 -*-
#!/usr/bin/env python
# @Project : bioinformatics
# @File : bwt2.py
# @Author : ycy
# @Time : 2023/6/3 22:25
# @Software : PyCharm Professional
"""

from time import time


def read_fasta(file_path):
    """
    Read the sequence from a FASTA file.
    """
    sequence = ''
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue  # Skip header lines starting with ">"
            sequence += line.strip()  # Append the current line to the sequence string, removing any leading/trailing whitespace.

    return sequence


def F_rank(bwt):
    """
    Compute the number of occurrences and first-occurrence positions of each character in the first column of the BWT matrix.
    """
    total_num = {}  # Dictionary to store the total number of occurrences of each character.
    ranks = []  # List to store the rank of each character in the last column.
    for i in bwt:
        if i not in total_num:
            total_num[i] = 0  # If character i is not in the dictionary yet, initialize its count to 0.
        ranks.append(total_num[i])  # Add the rank of the current character to the list of ranks.
        total_num[i] += 1  # Increment the count of character i.
    first = {}  # Dictionary to store the range of positions where each character appears in the first column.
    total_c = 0  # Total count of characters seen so far.
    for c, count in sorted(total_num.items()):  # Iterate over the sorted characters in the first column.
        first[c] = (total_c,
                    total_c + count)  # Store the range of positions where the current character appears in the first column.
        total_c += count  # Update the starting position for the next character.
    F_range = first
    return F_range, ranks


def bwt_encode(seq):
    """
    Encode a sequence using the Burrows-Wheeler transform (BWT).
    """
    # Remove whitespace characters from the sequence.
    seq = "".join(char for char in seq if not char.isspace())

    n = len(seq)
    # Append a special character to the end of the sequence.
    seq += "$"
    # Construct the matrix of rotations.
    rotations = [seq[i:] + seq[:i] for i in range(n + 1)]
    # Sort the matrix of rotations lexicographically.
    rotations.sort()
    # Get the offsets of the original order of rotations in the sorted matrix.
    offsets = [n - r.index("$") for r in rotations]
    # Get the last column of the sorted matrix.
    bwt = "".join([r[-1] for r in rotations])
    return bwt, offsets


def bwt_decode(bwt):
    """
    Decode a sequence from its Burrows-Wheeler transformed (BWT) representation.
    """
    str_F, ranks = F_rank(bwt)
    result = ""
    start_row = 0

    while bwt[start_row] != "$":
        c = bwt[start_row]  # Current character.
        result = c + result
        start_row = str_F[c][0] + ranks[start_row]
    return result


def bwt_search(bwt, query, offset):
    """
    Find the positions of a query sequence in the original sequence represented by the given Burrows-Wheeler transform (BWT).
    """
    # Remove whitespace characters from the query.
    query = "".join(char for char in query if not char.isspace())

    F_range, ranks = F_rank(bwt)
    dic = {}
    next = 0  # Used for indexing during search after the first character of the pattern.
    for i in bwt:
        dic[i] = dic.get(i, 0) + 1
    query = query[::-1]  # Reverse the query sequence.
    char = query[0]  # Get the first character of the reversed query.
    init_loc = F_range[char][0]  # Find the position of the first occurrence of char in the first column.
    for num in range(init_loc, init_loc + dic[char]):  # Iterate over all positions where the first character appears.
        count = 2  # Used to check if the search is successful.
        if len(query) == 1:
            start = offset[num]
            print(f'seq found in {start}--{start + len(query) - 1}')
            continue
        if bwt[num] == query[1]:  # Start matching.
            next = F_range[bwt[num]][0] + ranks[num]
            if len(query) == 2:
                start = offset[next]
                print(f'seq found in {start}--{start + len(query) - 1}')
                continue
        for i in range(2, len(query)):
            if bwt[next] == query[i]:
                next = F_range[bwt[next]][0] + ranks[next]
                count += 1
            if count == len(query):
                start = offset[next]
                print(f'seq found in {start}--{start + len(query) - 1}')


def main():
    # seq = "ACGTACGTGACG"
    seq = '''CTAAATGTTTACCCCATAGATGTGAAACAATGATTCTTCATATATTAACATATTTTTTGACTTATACTTTTC
    TTCATCTAGTAAGGCGTTAATTTTTTCCGGATCTGTCGTTTTTATTGATAAAAGAGAAGAGTCTGGACTGTA
    ATTTTTAAATAATAAGATATTTATTAATATCCAATTATTCGTTTGGCTCGCTATTTCCATGCTCTCTTCGAA
    AGCATCAGCTCCTAAATCTATACAAAGGAATAAGTTACCTTCACAAAAATTCATTACCGAGGTAATCATTGC
    CCGATTAATGTCAGCCCCCAACATAAAACAATAATATATAGTTGTATAATTACAATCATACATACAGGCCAA
    CTGCATCATTTCATCAATGTCTATATTTGTCTTCTCTTTGTTATAAATTTCATGAAGGTCAAAGACGTTGTT
    ATAAGCAACCCCACATATTAACCGCCAATCTTTAAAATGACTATATCGTTGATAAAAATATTGGATGGCTTC
    AGTAAGCTTATATAGTATCGCCATACTATACCAATACCTAGTTAGCATTTCGTTGAATGAAATATTATCCAA
    TGTAAAGTTAATTGATAATGTATCTAGTTCACCAAAAATTCTTAATTTCAGTTGAGCATTATTTAGGAAAAG
    GGGATTATCAGATAATAATTCATGGCATAGAATAATATTACTGCTAGTTTTAACATACTGTACATTATAAAA
    TATTTCTAAAATTTTATTTTCACTCAAAGCTTTCCTCGCACCTAACTTTTGGCATAGGTCCTGGTGCAC'''
    bwt, offset = bwt_encode(seq)
    print("seq:", seq)
    print("BWT:", bwt)
    orig_seq = bwt_decode(bwt)
    print("Original sequence:", orig_seq)
    query = 'TTTAAAATGACTATA'
    print("query:", query)
    t0 = time()
    bwt_search(bwt, query, offset)
    print("time:", time() - t0)


if __name__ == "__main__":
    main()
