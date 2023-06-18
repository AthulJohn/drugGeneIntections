from math import sqrt
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq


def smith_waterman_score(aminoseqs):
    # Define the amino acid sequences to be aligned
    score=[]
    for i in range(len(aminoseqs)):
        seq1 = Seq(aminoseqs[i])
        score.append([])
        for j in range(len(aminoseqs)):
            seq2 = Seq(aminoseqs[j])

            aligner = PairwiseAligner()
            aligner.mode = 'local'  # Set alignment mode to local (Smith-Waterman)
            aligner.match_score = 1  # Set match score to 1
            aligner.mismatch_score = -1  # Set mismatch score to -1
            aligner.open_gap_score = -1  # Set open gap penalty to -1
            aligner.extend_gap_score = -0.5 
            # Perform Smith-Waterman alignment
            alignments = aligner.align(seq1, seq2)

            # Extract the alignment score from the first alignment in the list
            alignment = alignments[0]
            score[i].append( alignment.score)

    normalised_score=[]
    for i in range(len(aminoseqs)):
        normalised_score.append([])
        for j in range(len(aminoseqs)):
            normalised_score[i].append(score[i][j]/(sqrt(score[i][i])*sqrt(score[i][j])))
            
    return normalised_score
