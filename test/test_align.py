# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    # Initializing the NeedlemanWunsch class
    nw = NeedlemanWunsch(
        sub_matrix_file="./substitution_matrices/BLOSUM62.mat",
        gap_open=-10,
        gap_extend=-1
    )
    # Performing the alignment
    _ = nw.align(seq1, seq2)

    # Drop '-' from seqA_align and seqB_align and assert that the sequences are equal to seq3 and seq4
    assert nw.seqA_align.replace("-", "") == seq1
    assert nw.seqB_align.replace("-", "") == seq2

    # Asserting that the alignment score is correct
    assert nw.alignment_score == 4

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    # Initializing the NeedlemanWunsch class
    nw = NeedlemanWunsch(
        sub_matrix_file="./substitution_matrices/BLOSUM62.mat",
        gap_open=-10,
        gap_extend=-1
    )
    # Performing the alignment
    nw.align(seq3, seq4)

    # Drop '-' from seqA_align and seqB_align and assert that the sequences are equal to seq3 and seq4
    assert nw.seqA_align.replace("-", "") == seq3
    assert nw.seqB_align.replace("-", "") == seq4

    # Asserting that the backtrace is correct
    assert nw.seqA_align == "MAVHQLIRRP"
    assert nw.seqB_align == "M---QLIRHP"
    assert nw.alignment_score == 17




