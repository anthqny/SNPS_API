def detectar_snps(ref, paciente):
    if len(ref) != len(paciente):
        raise ValueError("Las secuencias deben tener la misma longitud")

    snps = []

    for i in range(len(ref)):
        if ref[i] != paciente[i]:
            snps.append({
                "posicion": i + 1, 
                "referencia": ref[i],
                "mutacion": paciente[i]
            })

    return snps