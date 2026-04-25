import os

def read_fasta(file_name):
    if not os.path.exists(file_name):
        raise FileNotFoundError(f"No se encontró el archivo: {file_name}")

    if not file_name.endswith(".fasta"):
        raise ValueError("El archivo debe tener extensión .fasta")

    sequences = {}
    defline = None
    secuencia = []

    encontrado_header = False

    try:
        with open(file_name, "r") as f:
            for line in f:
                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):
                    encontrado_header = True

                    if defline:
                        sequences[defline] = "".join(secuencia)

                    defline = line[1:].strip()
                    secuencia = []

                else:
                    if defline is None:
                        raise ValueError("Formato FASTA inválido: secuencia sin encabezado")

                    if not all(base in "ATCG" for base in line.upper()):
                        raise ValueError(f"Caracter inválido en la secuencia: {line}")

                    secuencia.append(line.upper())

            if defline:
                sequences[defline] = "".join(secuencia)

    except Exception as e:
        raise RuntimeError(f"Error leyendo FASTA: {e}")

    if not encontrado_header:
        raise ValueError("El archivo no contiene encabezados FASTA (>)")

    if not sequences:
        raise ValueError("No se encontraron secuencias en el archivo")

    return sequences