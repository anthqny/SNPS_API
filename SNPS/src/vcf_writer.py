from datetime import datetime


def escribir_vcf(snps, output_path):
    with open(output_path, "w") as f:

        # -------------------------
        # HEADERS
        # -------------------------
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
        f.write("##source=MiVariantCaller_v3.0\n")

        f.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">\n')
        f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
        f.write('##INFO=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        # -------------------------
        # DATA
        # -------------------------
        for snp in snps:
            chrom = "chr1"
            pos = snp["posicion"]
            ref = snp["referencia"]
            alt = snp["mutacion"]

            qual = snp.get("qual", ".")

            info = (
                f"AF={snp['af']};"
                f"DP={snp['depth']};"
                f"GT={snp['genotipo']}"
            )

            line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\tPASS\t{info}\n"
            f.write(line)