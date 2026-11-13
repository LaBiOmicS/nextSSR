import random, sys
random.seed(int(sys.argv[1]))
motifs = ["A","C","AG","AT","GC","AAG","ACT","AAAT","AGCT","AAAAG","ACGTC","AAGAAC","ATCGGA","AC","CA","T","AAT"]
out = []
for r in range(int(sys.argv[2])):
    parts = []
    for _ in range(random.randint(0, 40)):
        parts.append("".join(random.choice("ACGT") for _ in range(random.choice([0,1,2,3,5,10,50,120,300]))))
        m = random.choice(motifs)
        n = random.randint(2, 25)
        rep = m * n
        if random.random() < 0.2: rep = rep[:-1]
        if random.random() < 0.15:
            i = random.randrange(len(rep)); rep = rep[:i] + random.choice("ACGTN") + rep[i+1:]
        if random.random() < 0.2: rep = rep.lower()
        if random.random() < 0.05: rep = "A"*random.randint(8,20) + rep
        parts.append(rep)
    seq = "".join(parts)
    if random.random() < 0.1: seq = seq.replace("C", "N", 3)
    hdr = f"seq{r}" + (f" sample desc {r}\tx" if r % 3 == 0 else "")
    eol = "\r\n" if r % 7 == 0 else "\n"
    w = random.choice([60, 70, 80, 10**9])
    lines = [seq[i:i+w] for i in range(0, len(seq), w)] or [""]
    out.append(">" + hdr + eol + eol.join(lines) + eol)
sys.stdout.write("".join(out))
