# quality.py

def phred_score(char):
    return ord(char) - 33


def phred_to_error_prob(q):
    return 10 ** (-q / 10)


def filtrar_calidad(seq, qual, threshold=20):
    nueva = []

    for base, q in zip(seq, qual):
        if phred_score(q) >= threshold:
            nueva.append(base)
        else:
            nueva.append("N")

    return "".join(nueva)


def calidad_promedio(qual):
    if not qual:
        return 0

    return sum(phred_score(q) for q in qual) / len(qual)


def porcentaje_buenas(qual, threshold=20):
    if not qual:
        return 0

    buenas = sum(1 for q in qual if phred_score(q) >= threshold)
    return buenas / len(qual)