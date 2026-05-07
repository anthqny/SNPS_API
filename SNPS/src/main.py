import os

from parser import read_fasta, validar_fasta_para_snp
from fastq_parser import read_fastq
from quality import filtrar_calidad
from aligner import alinear
from variant_counter import contar_variantes
from snp_detector import detectar_snps_con_af
from utils import guardar_json
from vcf_writer import escribir_vcf


def formatear_snps(nombre_archivo, snps):
    print(f"\nArchivo: {nombre_archivo}")

    if not snps:
        print("  No se encontraron SNPs")
        return

    for snp in snps:
        print(
            f"  🔹 Pos {snp['posicion']}: "
            f"{snp['referencia']}→{snp['mutacion']} "
            f"| AF={snp['af']} | DP={snp['depth']} "
            f"| QUAL={snp['qual']} | GT={snp['genotipo']}"
        )


def procesar_fastq(ruta_archivo, ref_global, nombre_archivo):
    reads = read_fastq(ruta_archivo)

    if not reads:
        raise ValueError("FASTQ vacío")

    reads_alineados = []
    ref_aligned_global = None

    for read in reads:
        seq_filtrada = filtrar_calidad(read["seq"], read["qual"])

        if not seq_filtrada or set(seq_filtrada) == {"N"}:
            continue

        ref_aligned, read_aligned = alinear(ref_global, seq_filtrada)

        ref_aligned_global = ref_aligned
        reads_alineados.append({
    "seq": read_aligned,
    "qual": read["qual"][:len(read_aligned)]
})

    if not reads_alineados:
        print(f"[WARN] {nombre_archivo}: sin lecturas útiles")
        return []

    conteo = contar_variantes(ref_aligned_global, reads_alineados)

    snps = detectar_snps_con_af(ref_aligned_global, conteo)

    return snps


def main():
    base_dir = os.path.dirname(os.path.dirname(__file__))

    input_dir = os.path.join(base_dir, "data", "input")
    output_dir = os.path.join(base_dir, "data", "output")

    os.makedirs(output_dir, exist_ok=True)

    archivos = os.listdir(input_dir)

    fasta_files = [f for f in archivos if f.endswith((".fasta", ".fa"))]
    fastq_files = [f for f in archivos if f.endswith(".fastq")]

    if not fasta_files:
        print("Necesitas un FASTA de referencia")
        return

    
    ref_path = os.path.join(input_dir, fasta_files[0])

    datos = read_fasta(ref_path)
    ref_global, _ = validar_fasta_para_snp(datos)

    print(f"\Procesando {len(fastq_files)} FASTQ...\n")

    resumen = []

    for archivo in fastq_files:
        ruta = os.path.join(input_dir, archivo)

        try:
            snps = procesar_fastq(ruta, ref_global, archivo)

            formatear_snps(archivo, snps)

            json_out = os.path.join(output_dir, archivo + ".json")
            vcf_out = os.path.join(output_dir, archivo + ".vcf")

            guardar_json(snps, json_out)
            escribir_vcf(snps, vcf_out)

            resumen.append((archivo, len(snps)))

        except Exception as e:
            print(f"Error en {archivo}: {e}")

    print("\nResumen final:")
    for a, n in resumen:
        print(f"  {a}: {n} SNP(s)")

    print("\n Pipeline terminado\n")


if __name__ == "__main__":
    main()