import os
from parser import read_fasta
from snp_detector import detectar_snps
from utils import guardar_json


def formatear_snps(nombre_archivo, snps):
    print(f"\n Archivo: {nombre_archivo}")

    if not snps:
        print("  ✔ No se encontraron SNPs")
        return

    for snp in snps:
        print(f"  🔹 Posición {snp['posicion']}: {snp['referencia']} → {snp['mutacion']}")


def main():
    try:
        base_dir = os.path.dirname(os.path.dirname(__file__))

        input_dir = os.path.join(base_dir, "data", "input")
        output_dir = os.path.join(base_dir, "data", "output")

        archivos = os.listdir(input_dir)

        fasta_files = [f for f in archivos if f.endswith(".fasta")]
        otros_archivos = [f for f in archivos if not f.endswith(".fasta")]

        if otros_archivos:
            print("\n⚠ Archivos ignorados (no FASTA):")
            for f in otros_archivos:
                print(f"   - {f}")

        if not fasta_files:
            print("\n No hay archivos FASTA para procesar")
            return

        print(f"\n Procesando {len(fasta_files)} archivos...\n")

        resumen_total = []

        for archivo in fasta_files:
            ruta_fasta = os.path.join(input_dir, archivo)

            nombre_salida = archivo.replace(".fasta", "_resultado.json")
            ruta_output = os.path.join(output_dir, nombre_salida)

            try:
                datos = read_fasta(ruta_fasta)

                ref = datos.get("ref")
                paciente = datos.get("paciente")

                if not ref or not paciente:
                    raise ValueError("Faltan secuencias 'ref' y 'paciente'")

                snps = detectar_snps(ref, paciente)

                
                formatear_snps(archivo, snps)

            
                guardar_json({
                    "archivo": archivo,
                    "total_snps": len(snps),
                    "snps": snps
                }, ruta_output)

                resumen_total.append({
                    "archivo": archivo,
                    "total_snps": len(snps)
                })

            except Exception as e:
                print(f"\n Error en {archivo}: {e}")


        print("\n Resumen final:")
        for r in resumen_total:
            print(f"  {r['archivo']}: {r['total_snps']} SNP(s)")

        print("\n Proceso terminado correctamente.\n")

    except Exception as e:
        print(f"\n Error general: {e}")


if __name__ == "__main__":
    main()