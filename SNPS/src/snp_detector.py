import math


def calcular_qual(af, avg_q):

    if af >= 1:
        af_component = 60
    else:
        error_prob = 1 - af

        if error_prob <= 0:
            af_component = 60
        else:
            af_component = -10 * math.log10(error_prob)

    # combinar AF + calidad promedio
    qual = af_component * (avg_q / 40)

    return round(min(qual, 60), 2)


def clasificar_genotipo(af):
    if af > 0.8:
        return "HOM"
    elif af > 0.3:
        return "HET"
    return "LOW"


def detectar_snps_con_af(
    ref,
    conteo,
    threshold=0.2,
    min_dp=5,
    min_qual=10 #Aca se cambia para ser mas estricto 10 o 20 super estricto; 5 mas sensible !Sensitivity vs Specificity
):

    snps = []

    for i, ref_base in enumerate(ref[:len(conteo)]):

        total = sum(conteo[i]["bases"].values())

        if total < min_dp:
            continue

        for base, count in conteo[i]["bases"].items():

            if base == ref_base:
                continue

            af = count / total

            if af < threshold:
                continue

            avg_q = (
                conteo[i]["qual_sum"][base] / count
            )

            qual = calcular_qual(af, avg_q)

            if qual < min_qual:
                continue

            gt = clasificar_genotipo(af)

            snps.append({
                "posicion": i + 1,
                "referencia": ref_base,
                "mutacion": base,
                "depth": total,
                "alt_count": count,
                "af": round(af, 3),
                "avg_q": round(avg_q, 2),
                "qual": qual,
                "genotipo": gt
            })

    return snps