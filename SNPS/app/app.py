from fastapi import FastAPI, Request, UploadFile, File
from fastapi.templating import Jinja2Templates
from fastapi.responses import HTMLResponse
from Bio import SeqIO
import math
import io

app = FastAPI()

templates = Jinja2Templates(
    directory="app/templates"
)

MIN_PHRED = 10


def calcular_qual(af, dp):

    if af == 1:
        return 60.0

    error_prob = max(1e-6, 1 - af)

    qual = -10 * math.log10(error_prob) * dp / 10

    return round(qual, 2)


@app.get("/", response_class=HTMLResponse)
async def home(request: Request):

    return templates.TemplateResponse(
        request=request,
        name="index.html",
        context={
            "resultados": None
        }
    )


@app.post("/upload", response_class=HTMLResponse)
async def upload_files(
    request: Request,
    fasta: UploadFile = File(...),
    fastq: UploadFile = File(...)
):

    resultados = []

    # =========================
    # LEER FASTA
    # =========================

    fasta_content = await fasta.read()

    fasta_io = io.StringIO(
        fasta_content.decode("utf-8")
    )

    fasta_records = list(
        SeqIO.parse(fasta_io, "fasta")
    )

    referencia = str(
        fasta_records[0].seq
    )

    # =========================
    # LEER FASTQ
    # =========================

    fastq_content = await fastq.read()

    fastq_io = io.StringIO(
        fastq_content.decode("utf-8")
    )

    reads = list(
        SeqIO.parse(fastq_io, "fastq")
    )

    conteo = {}

    total_reads = 0

    # =========================
    # COMPARAR
    # =========================

    for read in reads:

        seq = str(read.seq)

        quals = read.letter_annotations[
            "phred_quality"
        ]

        total_reads += 1

        for i in range(
            min(len(seq), len(referencia))
        ):

            if quals[i] < MIN_PHRED:
                continue

            ref_base = referencia[i]

            alt_base = seq[i]

            if ref_base == alt_base:
                continue

            key = (
                i + 1,
                ref_base,
                alt_base
            )

            conteo[key] = (
                conteo.get(key, 0) + 1
            )

    # =========================
    # RESULTADOS
    # =========================

    for (pos, ref, alt), count in conteo.items():

        af = round(
            count / total_reads,
            2
        )

        gt = (
            "HOM"
            if af >= 0.8
            else "HET"
        )

        qual = calcular_qual(
            af,
            total_reads
        )

        resultados.append({

            "pos": pos,
            "ref": ref,
            "alt": alt,
            "af": af,
            "dp": total_reads,
            "qual": qual,
            "gt": gt

        })

    return templates.TemplateResponse(
        request=request,
        name="index.html",
        context={
            "resultados": resultados
        }
    )