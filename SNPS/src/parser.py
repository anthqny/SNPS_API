import os

# Bases válidas según estándar IUPAC
IUPAC_BASES = set("ATCGRYSWKMBDHVN")


def read_fasta(file_name):
    """
    Lee un archivo FASTA y devuelve un diccionario:
    {
        "header1": "SECUENCIA",
        "header2": "SECUENCIA"
    }
    """

    if not os.path.exists(file_name):
        raise FileNotFoundError(f"No se encontró el archivo: {file_name}")

    if not file_name.lower().endswith((".fasta", ".fa", ".fna")):
        raise ValueError("El archivo debe ser formato FASTA (.fasta, .fa, .fna)")

    sequences = {}
    defline = None
    secuencia = []

    try:
        with open(file_name, "r") as f:
            for line_num, line in enumerate(f, start=1):
                line = line.strip()

                # Ignorar líneas vacías
                if not line:
                    continue

                # Header
                if line.startswith(">"):
                    if defline:
                        sequences[defline] = "".join(secuencia)

                    defline = line[1:].strip()

                    if not defline:
                        raise ValueError(f"Header vacío en línea {line_num}")

                    secuencia = []

                else:
                    if defline is None:
                        raise ValueError(
                            f"Formato FASTA inválido: secuencia sin encabezado (línea {line_num})"
                        )

                    line = line.upper()

                    # Validar caracteres IUPAC
                    invalid_chars = set(line) - IUPAC_BASES
                    if invalid_chars:
                        raise ValueError(
                            f"Caracteres inválidos {invalid_chars} en línea {line_num}"
                        )

                    secuencia.append(line)

            # Guardar última secuencia
            if defline:
                sequences[defline] = "".join(secuencia)

    except Exception as e:
        raise RuntimeError(f"Error leyendo FASTA: {e}")

    # Validaciones finales
    if not sequences:
        raise ValueError("No se encontraron secuencias en el archivo FASTA")

    # Debug opcional
    print(f"[INFO] Secuencias cargadas: {list(sequences.keys())}")

    return sequences


def validar_fasta_para_snp(sequences):
    """
    Valida que existan al menos dos secuencias
    y que sean comparables
    """

    if len(sequences) < 2:
        raise ValueError("Se requieren al menos 2 secuencias para comparar SNPs")

    nombres = list(sequences.keys())

    ref = sequences[nombres[0]]
    paciente = sequences[nombres[1]]

    if len(ref) == 0 or len(paciente) == 0:
        raise ValueError("Las secuencias no pueden estar vacías")

    return ref, paciente