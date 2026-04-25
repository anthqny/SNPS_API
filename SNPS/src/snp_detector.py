from quality import phred_score, phred_to_error_prob


def detectar_snps(ref, paciente, qual=None):
    snps = []

    for i, (r, p) in enumerate(zip(ref, paciente)):
        if r == "-" or p == "-":
            continue

        if r == "N" or p == "N":
            continue

        if r != p:
            snp = {
                "posicion": i + 1,
                "referencia": r,
                "mutacion": p,
                "tipo": "SNP"
            }

            # 🔥 NUEVO: calcular confianza si hay calidad
            if qual:
                q_score = phred_score(qual[i])
                error_prob = phred_to_error_prob(q_score)
                confianza = 1 - error_prob

                snp["phred"] = q_score
                snp["confianza"] = round(confianza, 6)

            snps.append(snp)

    return snps