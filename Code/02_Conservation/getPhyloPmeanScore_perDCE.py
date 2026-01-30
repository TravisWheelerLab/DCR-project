import pyBigWig
import pandas as pd

bw = pyBigWig.open("./data/hg38.phyloP100way.bw")
per_exon = []
for line in open("output/02_Conservation/DCE_CDS_intersect_filtered.bed"): # here there will be duplicate exons from DCR that share the same exons. It is not important because we only report per-DCR score and find outliers.
    cols = line.strip().split()
    mean = bw.stats(cols[0], int(cols[1]), int(cols[2])) # get the value field
    exon = cols[3]
    eLen = int(cols[2])-int(cols[1])
    dce = cols[4]
    #print(mean[0])
    if mean[0] is None:
        print("WARNING: the PP score of exon " + exon + " of DCE: " + dce + " cannot be calculated")
    else:
        weighted_mean_PP = (mean[0]) * eLen
        output = [exon, mean[0], dce, eLen, weighted_mean_PP]
        per_exon.append(output)
bw.close()

per_exon_table = pd.DataFrame(per_exon, columns=['exon', 'mean_PP', 'DCEindex', 'length', 'weighted_mean_PP'])
DCEuniq = pd.unique(per_exon_table['DCEindex'])
per_dce = []
for DCE in DCEuniq:
    DCEtable = per_exon_table[per_exon_table['DCEindex'] == DCE]    
    DCEpp = sum(DCEtable['weighted_mean_PP'])/sum(DCEtable['length'])
    per_dce.append([DCE, DCEpp])

final_table = pd.DataFrame(per_dce , columns=['DCEindex', 'mean_PP'])
final_table.to_csv('output/02_Conservation/meanPP_per_DCE.tsv', sep = "\t", index = False)