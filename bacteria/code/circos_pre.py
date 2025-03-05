from Bio import SeqIO
recs = SeqIO.parse("./outdir/Prokka/sample/sample.gbk", "genbank")
# 使用 SeqIO.parse() 来处理多个记录
forward = []
reverse = []
total = 0
for rec in recs:
    length = len(rec)
    for feature in rec.features:
        if feature.type != 'CDS':
           continue
        if feature.location.strand == 1:
            forward.append((feature.qualifiers['locus_tag'][0], int(feature.location.start) + total, int(feature.location.end) + total))
        elif feature.location.strand == -1:
            reverse.append((feature.qualifiers['locus_tag'][0], int(feature.location.start) + total, int(feature.location.end) + total))
    total += length
def get_forward_strand(strands):
    """
    Create forward strands
    """
    with open("sample/forward.txt", "w") as f:
        for strand in strands:
            f.writelines("chr1\t%d\t%d\t%s\n" % (strand[1], strand[2], strand[0]))

def get_reverse_strand(strands):
    """
    Create reverse strands
    """
    with open("sample/reverse.txt", "w") as f:
        for strand in strands:
            f.writelines("chr1\t%d\t%d\t%s\n" % (strand[1], strand[2], strand[0]))
def calc_GC_content(seq):
    """
    Return GC content of input sequence
    """
    gc = sum(seq.count(x) for x in ['G','C','g','c'])
    gc_content = gc/float(len(seq))
    return round(gc_content, 4)

def calc_GC_skew(seq):
    """
    Reuturn GC skew of input sequence
    """
    g = seq.count('G') + seq.count('g')
    c = seq.count('C') + seq.count('c')
    try:
        skew = (g - c)/float(g + c)
    except ZeroDivisionError:
        skew = 0
    return round(skew, 4)

window = 10000
step = 5000
seqs = SeqIO.parse("./outdir/Prokka/sample/sample.fna", "fasta")
seq_total = 0
with open("sample/gc_skew.txt", "w") as f1, open("sample/gc_content.txt", "w") as f2:
    for seq in seqs:
        seq_1 = seq.seq
        seq_len = len(seq_1)
        for i in range(0, len(seq_1), step):
            subseq = seq_1[i:i+window]
            gc_content = calc_GC_content(subseq)
            gc_skew = calc_GC_skew(subseq)
            i += seq_total
            start = i + 1 if i + 1 <= len(seq_1)+seq_total else i
            end = i + step if i + step <= len(seq_1)+seq_total else len(seq_1)+seq_total
            f1.writelines("chr1\t%d\t%d\t%s\n" % (start, end, gc_skew))
            f2.writelines("chr1\t%d\t%d\t%s\n" % (start, end, gc_content))
        seq_total += seq_len
get_forward_strand(forward)
get_reverse_strand(reverse)