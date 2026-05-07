import math


# -------------------------
# CALCULAR QUAL (Phred)
# -------------------------
def calcular_qual(af):
    """
    Convierte frecuencia de variante en score Phred
    """
    if af >= 1:
        return 60  # cap estándar

    error_prob = 1 - af

    if error_prob <= 0:
        return 60

    return round(-10 * math.log10(error_prob), 2)


# -------------------------
# CLASIFICAR GENOTIPO
# -------------------------
def clasificar_genotipo(af):
    """
    Clasificación simple tipo diploide
    """
    if af > 0.8:
        return "HOM"  # homocigoto
    elif af > 0.3:
        return "HET"  # heterocigoto
    else:
        return "LOW"  # baja frecuencia


# -------------------------
# DETECTOR PRINCIPAL
# -------------------------
def detectar_snps_con_af(
    ref,
    conteo,
    threshold=0.2,
    min_dp=5
):
    """
    SNP caller basado en:
    - AF mínimo
    - profundidad mínima
    - cálculo de QUAL
    """

    snps = []

    for i, ref_base in enumerate(ref[:len(conteo)]):

        total = sum(conteo[i].values())

        # ❌ filtro de profundidad
        if total < min_dp:
            continue

        for base, count in conteo[i].items():

            # ❌ ignorar referencia
            if base == ref_base:
                continue

            af = count / total

            # ❌ filtro AF
            if af < threshold:
                continue

            qual = calcular_qual(af)
            gt = clasificar_genotipo(af)

            snps.append({
                "posicion": i + 1,
                "referencia": ref_base,
                "mutacion": base,
                "depth": total,
                "alt_count": count,
                "af": round(af, 3),
                "qual": qual,
                "genotipo": gt
            })

    return snps