from Bio.Align import PairwiseAligner

def alinear(ref, paciente):
    aligner = PairwiseAligner()

    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1

    alignment = aligner.align(ref, paciente)[0]

    return str(alignment.target), str(alignment.query)