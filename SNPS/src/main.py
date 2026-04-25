import os

from parser import read_fasta, validar_fasta_para_snp
from fastq_parser import read_fastq
from quality import filtrar_calidad, calidad_promedio, porcentaje_buenas
from aligner import alinear
from snp_detector import detectar_snps
from utils import guardar_json


def formatear_snps(nombre_archivo, snps):
    print(f"\n📄 Archivo: {nombre_archivo}")

    if not snps:
        print("  ✔ No se encontraron SNPs")
        return

    for snp in snps:
        linea = (
            f"  🔹 Posición {snp['posicion']}: "
            f"{snp['referencia']} → {snp['mutacion']}"
        )

        if "confianza" in snp:
            linea += f" | Q={snp['phred']} | conf={snp['confianza']}"

        print(linea)


# -------------------------
# FASTA
# -------------------------
def procesar_fasta(ruta_archivo):
    datos = read_fasta(ruta_archivo)
    ref, paciente = validar_fasta_para_snp(datos)

    ref_aligned, paciente_aligned = alinear(ref, paciente)

    snps = detectar_snps(ref_aligned, paciente_aligned)

    return snps


# -------------------------
# FASTQ
# -------------------------
def procesar_fastq(ruta_archivo, ref_global, nombre_archivo):
    reads = read_fastq(ruta_archivo)

    if not reads:
        raise ValueError("El archivo FASTQ no contiene lecturas")

    # MVP: usar solo el primer read
    read = reads[0]

    seq_filtrada = filtrar_calidad(read["seq"], read["qual"])
    calidad = calidad_promedio(read["qual"])
    porcentaje = porcentaje_buenas(read["qual"])

    ref_aligned, paciente_aligned = alinear(ref_global, seq_filtrada)

    snps = detectar_snps(ref_aligned, paciente_aligned, read["qual"])

    # ⚠ warning bien contextualizado
    if calidad < 20:
        print(f"[WARN] {nombre_archivo}: Baja calidad ({calidad:.2f})")

    return snps, calidad, porcentaje


# -------------------------
# MAIN
# -------------------------
def main():
    try:
        base_dir = os.path.dirname(os.path.dirname(__file__))

        input_dir = os.path.join(base_dir, "data", "input")
        output_dir = os.path.join(base_dir, "data", "output")

        archivos = os.listdir(input_dir)

        fasta_files = [f for f in archivos if f.lower().endswith((".fasta", ".fa", ".fna"))]
        fastq_files = [f for f in archivos if f.lower().endswith(".fastq")]
        otros_archivos = [
            f for f in archivos
            if f not in fasta_files and f not in fastq_files
        ]

        if otros_archivos:
            print("\n⚠ Archivos ignorados:")
            for f in otros_archivos:
                print(f"   - {f}")

        if not fasta_files and not fastq_files:
            print("\n❌ No hay archivos FASTA o FASTQ para procesar")
            return

        print(f"\n🚀 Procesando {len(fasta_files)} FASTA y {len(fastq_files)} FASTQ...\n")

        resumen_total = []

        # -------------------------
        # FASTA
        # -------------------------
        for archivo in fasta_files:
            ruta = os.path.join(input_dir, archivo)
            nombre_salida = archivo + "_resultado.json"
            ruta_output = os.path.join(output_dir, nombre_salida)

            try:
                snps = procesar_fasta(ruta)

                formatear_snps(archivo, snps)

                guardar_json({
                    "archivo": archivo,
                    "tipo": "FASTA",
                    "total_snps": len(snps),
                    "snps": snps
                }, ruta_output)

                resumen_total.append({
                    "archivo": archivo,
                    "tipo": "FASTA",
                    "total_snps": len(snps)
                })

            except Exception as e:
                print(f"\n❌ Error en {archivo}: {e}")

        # -------------------------
        # FASTQ
        # -------------------------
        if fasta_files:
            ref_path = os.path.join(input_dir, fasta_files[0])
            datos_ref = read_fasta(ref_path)
            ref_global, _ = validar_fasta_para_snp(datos_ref)

            for archivo in fastq_files:
                ruta = os.path.join(input_dir, archivo)
                nombre_salida = archivo + "_resultado.json"
                ruta_output = os.path.join(output_dir, nombre_salida)

                try:
                    snps, calidad, porcentaje = procesar_fastq(
                        ruta, ref_global, archivo
                    )

                    formatear_snps(archivo, snps)

                    guardar_json({
                        "archivo": archivo,
                        "tipo": "FASTQ",
                        "total_snps": len(snps),
                        "calidad_promedio": round(calidad, 2),
                        "porcentaje_bases_buenas": round(porcentaje, 2),
                        "snps": snps
                    }, ruta_output)

                    resumen_total.append({
                        "archivo": archivo,
                        "tipo": "FASTQ",
                        "total_snps": len(snps),
                        "calidad_promedio": round(calidad, 2),
                        "porcentaje_bases_buenas": round(porcentaje, 2)
                    })

                except Exception as e:
                    print(f"\n❌ Error en {archivo}: {e}")

        else:
            if fastq_files:
                print("\n⚠ No hay FASTA de referencia para procesar FASTQ")

        # -------------------------
        # RESUMEN FINAL
        # -------------------------
        print("\n📊 Resumen final:")
        for r in resumen_total:
            if r["tipo"] == "FASTQ":
                print(
                    f"  {r['archivo']} ({r['tipo']}): "
                    f"{r['total_snps']} SNP(s), "
                    f"calidad={r['calidad_promedio']}, "
                    f"bases_buenas={r['porcentaje_bases_buenas']}"
                )
            else:
                print(
                    f"  {r['archivo']} ({r['tipo']}): "
                    f"{r['total_snps']} SNP(s)"
                )

        print("\n✅ Proceso terminado correctamente.\n")

    except Exception as e:
        print(f"\n🔥 Error general: {e}")


if __name__ == "__main__":
    main()