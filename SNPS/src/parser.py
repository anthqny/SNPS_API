import os


def read_fasta(file_name):
    if not os.path.exists(file_name):
        raise FileNotFoundError(f"No se encontró el archivo: {file_name}")

    if not file_name.endswith((".fasta", ".fa")):
        raise ValueError("El archivo debe ser FASTA (.fasta o .fa)")

    sequences = {}
    defline = None
    secuencia = []

    try:
        with open(file_name, "r") as f:
            for line in f:
                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):
                    if defline:
                        sequences[defline] = "".join(secuencia)

                    defline = line[1:].strip()
                    secuencia = []

                else:
                    if defline is None:
                        raise ValueError("FASTA inválido: secuencia sin encabezado")

                    secuencia.append(line.upper())

            if defline:
                sequences[defline] = "".join(secuencia)

    except Exception as e:
        raise RuntimeError(f"Error leyendo FASTA: {e}")

    if not sequences:
        raise ValueError("No se encontraron secuencias en FASTA")

    print(f"[INFO] Secuencias cargadas: {list(sequences.keys())}")

    return sequences


# 🔥 NUEVA FUNCIÓN (CLAVE)
def validar_fasta_para_snp(datos):
    """
    Pipeline moderno:
    - solo usamos UNA referencia
    """

    if not datos:
        raise ValueError("FASTA vacío")

    # tomar la primera secuencia como referencia
    ref = list(datos.values())[0]

    return ref, None