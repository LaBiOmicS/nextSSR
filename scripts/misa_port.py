"""Line-by-line Python port of misa.pl v2.0 (non-GFF mode) used as a test oracle."""
import re, sys
from collections import defaultdict

def perl_substr(s, off, ln):
    n = len(s)
    if ln >= 0: return s[off:off+ln]
    end = n + ln
    return s[off:end] if end > off else ""

def main(fasta, ini):
    typrep = {}; amb = 0
    for line in open(ini):
        m = re.match(r'def\S*\s+(.*)', line, re.I)
        if m:
            nums = [int(x) for x in re.findall(r'(\d+)', m.group(1))]
            typrep = dict(zip(nums[0::2], nums[1::2]))
        m = re.match(r'int\S*\s+(\d+)', line, re.I)
        if m: amb = int(m.group(1))
    typ = sorted(typrep)
    out = open(fasta + ".misa.port", "w")
    out.write("ID\tSSR nr.\tSSR type\tSSR\tsize\tstart\tend\n")
    data = open(fasta, encoding="latin-1", newline="").read()
    max_repeats = 1; min_repeats = 1000
    count_motifs = defaultdict(int); motif_rep = defaultdict(lambda: defaultdict(int)); count_class = defaultdict(int)
    number_sequences = 0; size_sequences = 0; ssr_containing = defaultdict(int); ssr_in_compound = 0
    chunks = data.split(">")
    recs = [c + ">" for c in chunks[:-1]] + [chunks[-1]]
    for rec in recs:
        m = re.match(r'(.*?)\n(.*)', rec, re.S)
        if not m: continue
        id_, seq = m.group(1), m.group(2)
        seq = re.sub(r'[\d\s>]', '', seq)
        id_ = re.sub(r'^\s*', '', id_); id_ = re.sub(r'\s*$', '', id_); id_ = re.sub(r'\s', '_', id_)
        number_sequences += 1; size_sequences += len(seq)
        ssrs = []
        for k in typ:
            minreps = typrep[k] - 1
            if min_repeats > typrep[k]: min_repeats = typrep[k]
            for mt in re.finditer(r'(([acgt]{%d})\2{%d,})' % (k, minreps), seq, re.I):
                motif = mt.group(2).upper()
                red = any(len(motif) % j == 0 and re.search(r'([ACGT]{%d})\1{%d}' % (j, k // j - 1), motif) for j in range(k - 1, 0, -1))
                if red: continue
                ssr = mt.group(1); reps = len(ssr) // k; end = mt.end(); start = end - len(ssr) + 1
                ssrs.append((start, end, motif, reps))
                count_motifs[motif] += 1; motif_rep[motif][reps] += 1; count_class[k] += 1
                max_repeats = max(max_repeats, reps)
        if not ssrs: continue
        ssr_containing[len(ssrs)] += 1
        ssrs.sort(key=lambda x: x[0])
        one = lambda s: "(%s)%d" % (s[2], s[3])
        i = 0; cnt = 0; nr = len(ssrs)
        while i < nr:
            if i + 1 >= nr or ssrs[i+1][0] - ssrs[i][1] > amb + 1:
                cnt += 1; s = ssrs[i]
                out.write("%s\t%d\tp%d\t%s\t%d\t%d\t%d\n" % (id_, cnt, len(s[2]), one(s), s[1]-s[0]+1, s[0], s[1]))
                i += 1; continue
            cnt += 1; ssr_in_compound += 1
            a, b = ssrs[i], ssrs[i+1]
            if b[0] - a[1] == 1:
                ty = "c*"; pat = one(a) + one(b) + "*"
            else:
                ty = "c"; pat = one(a) + perl_substr(seq, a[1], b[0]-a[1]-1).lower() + one(b)
            start = a[0]; end = b[1]
            i += 1
            while i + 1 < nr and ssrs[i+1][0] - ssrs[i][1] <= amb + 1:
                ssr_in_compound += 1
                a, b = ssrs[i], ssrs[i+1]
                if b[0] - a[1] == 1:
                    pat += one(b) + "*"; ty = "c*"
                else:
                    pat += perl_substr(seq, a[1], b[0]-a[1]-1).lower() + one(b)
                end = b[1]; i += 1
            i += 1
            out.write("%s\t%d\t%s\t%s\t%d\t%d\t%d\n" % (id_, cnt, ty, pat, end-start+1, start, end))
    out.close()
    o = open(fasta + ".statistics.port", "w")
    o.write('Specifications\n==============\n\nSequence source file: "%s"\n\nDefinement of microsatellites (unit size / minimum number of repeats):\n' % fasta)
    o.write("".join("(%d/%d) " % (k, typrep[k]) for k in typ) + "\n")
    if amb > 0: o.write("\nMaximal number of bases interrupting 2 SSRs in a compound microsatellite:  %d\n" % amb)
    o.write("\n\n\n")
    total = sum(count_class.values()) if count_class else ""
    sc = sum(ssr_containing.values())
    o.write("RESULTS OF MICROSATELLITE SEARCH\n================================\n\n")
    o.write("Total number of sequences examined:              %d\n" % number_sequences)
    o.write("Total size of examined sequences (bp):           %d\n" % size_sequences)
    o.write("Total number of identified SSRs:                 %s\n" % total)
    o.write("Number of SSR containing sequences:              %d\n" % sc)
    o.write("Number of sequences containing more than 1 SSR:  %d\n" % (sc - ssr_containing.get(1, 0)))
    o.write("Number of SSRs present in compound formation:    %d\n\n\n" % ssr_in_compound)
    o.write("Distribution to different repeat type classes\n---------------------------------------------\n\nUnit size\tNumber of SSRs\n")
    for k in sorted(count_class): o.write("%d\t%d\n" % (k, count_class[k]))
    o.write("\n")
    o.write("Frequency of identified SSR motifs\n----------------------------------\n\nRepeats")
    for i in range(min_repeats, max_repeats + 1): o.write("\t%d" % i)
    o.write("\ttotal\n")
    key = lambda s: (len(s), s)
    for mo in sorted(count_motifs, key=key):
        o.write(mo)
        for j in range(min_repeats, max_repeats + 1):
            if j < typrep[len(mo)]: o.write("\t-"); continue
            o.write("\t%d" % motif_rep[mo][j] if motif_rep[mo].get(j) else "\t")
        o.write("\t%d\n" % count_motifs[mo])
    o.write("\n")
    o.write("Frequency of classified repeat types (considering sequence complementary)\n-------------------------------------------------------------------------\n\nRepeats")
    comp = str.maketrans("ACGT", "TGCA")
    groups = {}; seen = set()
    for mo in sorted(count_motifs, key=key):
        if mo in seen: continue
        act = acta = mo; rev = reva = mo.translate(comp)[::-1]; grp = set()
        for _ in range(len(mo)):
            if act in count_motifs: grp.add(act); seen.add(act)
            if rev in count_motifs: grp.add(rev); seen.add(rev)
            act = act[1:] + act[0]; rev = rev[1:] + rev[0]
            if act < acta: acta = act
            if rev < reva: reva = rev
        name = "%s/%s" % (acta, reva) if acta < reva else "%s/%s" % (reva, acta)
        g = groups.setdefault(name, defaultdict(int))
        for mm in grp:
            for kk in range(min_repeats, max_repeats + 1):
                if motif_rep[mm].get(kk): g["total"] += motif_rep[mm][kk]; g[kk] += motif_rep[mm][kk]
    for i in range(min_repeats, max_repeats + 1): o.write("\t%d" % i)
    o.write("\ttotal\n")
    for name in sorted(groups, key=key):
        t = (len(name) - 1) // 2; o.write(name)
        for j in range(min_repeats, max_repeats + 1):
            if j < typrep[t]: o.write("\t-"); continue
            o.write("\t%d" % groups[name][j] if groups[name].get(j) else "\t")
        o.write("\t%s\n" % groups[name].get("total", ""))

main(sys.argv[1], sys.argv[2])
