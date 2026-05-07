from collections import defaultdict

from quality import phred_score


def contar_variantes(ref, reads_alineados):
    
    if not reads_alineados:
        return []

    min_len = min(len(ref), *(len(r["seq"]) for r in reads_alineados))

    conteo = []

    for _ in range(min_len):
        conteo.append({
            "bases": defaultdict(int),
            "qual_sum": defaultdict(float)
        })

    for read in reads_alineados:

        seq = read["seq"]
        qual = read["qual"]

        for i in range(min_len):

            base = seq[i]

            if base in ["N", "-"]:
                continue

            q = phred_score(qual[i])

            conteo[i]["bases"][base] += 1
            conteo[i]["qual_sum"][base] += q

    return conteo