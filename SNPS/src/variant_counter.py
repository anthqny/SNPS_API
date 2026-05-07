from collections import defaultdict


def contar_variantes(ref, reads_alineados):
    """
    Cuenta bases ignorando ruido y manejando longitudes distintas
    """

    if not reads_alineados:
        return []

    min_len = min(len(ref), *(len(r) for r in reads_alineados))

    conteo = [defaultdict(int) for _ in range(min_len)]

    for read in reads_alineados:
        for i in range(min_len):
            base = read[i]

            # ❌ ignorar basura
            if base in ["-", "N"]:
                continue

            conteo[i][base] += 1

    return conteo