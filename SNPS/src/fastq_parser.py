def read_fastq(file_name):
    reads = []

    with open(file_name, "r") as f:
        while True:
            header = f.readline().strip()
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()

            if not qual:
                break

            reads.append({
                "header": header,
                "seq": seq.upper(),
                "qual": qual
            })

    return reads