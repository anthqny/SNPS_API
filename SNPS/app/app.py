from fastapi import FastAPI, Request, UploadFile, File
from fastapi.templating import Jinja2Templates
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles

from Bio import SeqIO

import math
import io

app = FastAPI()

app.mount(
    "/static",
    StaticFiles(directory="app/static"),
    name="static"
)

templates = Jinja2Templates(directory="app/templates")

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
        request,
        "index.html",
        {
            "request": request,
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

    fastq_content = await fastq.read()

    fastq_io = io.StringIO(
        fastq_content.decode("utf-8")
    )

    reads = list(
        SeqIO.parse(fastq_io, "fastq")
    )

    conteo = {}

    total_reads = len(reads)

    for read in reads:

        seq = str(read.seq)

        quals = read.letter_annotations["phred_quality"]

        for i in range(min(len(seq), len(referencia))):

            if quals[i] < MIN_PHRED:
                continue

            ref_base = referencia[i]

            alt_base = seq[i]

            if ref_base == alt_base:
                continue

            key = (i + 1, ref_base, alt_base)

            if key not in conteo:

                conteo[key] = {
                    "count": 0,
                    "read_seq": seq
                }

            conteo[key]["count"] += 1

    for (pos, ref, alt), data in conteo.items():

        count = data["count"]

        af = round(count / total_reads, 2)

        gt = "HOM" if af >= 0.8 else "HET"

        qual = calcular_qual(af, total_reads)

        read_seq = data["read_seq"]

        pointer = " " * (pos + 6) + "^"

        resultados.append({

            "pos": pos,
            "ref": ref,
            "alt": alt,
            "af": af,
            "dp": total_reads,
            "qual": qual,
            "gt": gt,

            "referencia": referencia,
            "read": read_seq,
            "pointer": pointer
        })

    return templates.TemplateResponse(
        request,
        "index.html",
        {
            "request": request,
            "resultados": resultados
        }
    )