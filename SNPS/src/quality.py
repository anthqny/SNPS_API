def phred_score(char):
    return ord(char) - 33


def prob_error(q):
    return 10 ** (-q / 10)


def filtrar_calidad(seq, qual, min_q=20):

    resultado = []

    for base, qchar in zip(seq, qual):

        q = phred_score(qchar)

        if q < min_q:
            resultado.append("N")
        else:
            resultado.append(base)

    return "".join(resultado)

def promedio_calidad(qual):
    scores = [phred_score(q) for q in qual]

    if not scores:
        return 0

    return sum(scores) / len(scores)