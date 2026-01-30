import pyBigWig
import sys

if len(sys.argv) < 2:
    print("Usage: python getPhyloPmeanScore.py <bigwig_file>")
    sys.exit(1)

bw = pyBigWig.open(sys.argv[1])

for line in open("output/02_Conservation/DCR_exons_DCR_genes_UNIQ.bed"):
    cols = line.strip().split()
    mean = bw.stats(cols[0], int(cols[1]), int(cols[2])) # get the value field
    exon = cols[3]
    eLen = int(cols[2])-int(cols[1])
    output = exon + "\t" + str(mean[0]) + "\tDCE_exon\t" + str(eLen)
    print(output)# Do something with the values...

for line in open("output/02_Conservation/nonDCR_exon_DCR_gene_UNIQ.bed"):
    cols = line.strip().split()
    try:
        mean = bw.stats(cols[0], int(cols[1]), int(cols[2])) # get the value field
        exon = cols[3]
        eLen = int(cols[2])-int(cols[1])
        output = exon + "\t" + str(mean[0]) + "\tnDCE_exon_DCE_gene\t" + str(eLen)
        print(output)# Do something with the values...
    except RuntimeError:
        pass

for line in open("output/02_Conservation/nonDCR_exon_nonDCR_gene_UNIQ.bed"):
    cols = line.strip().split()
    try:
        mean = bw.stats(cols[0], int(cols[1]), int(cols[2])) # get the value field
        exon = cols[3]
        eLen = int(cols[2])-int(cols[1])
        output = exon + "\t" + str(mean[0]) + "\tnDCE_exon_nDCE_gene\t" + str(eLen)
        print(output)# Do something with the values...
    except RuntimeError:
        pass
bw.close()