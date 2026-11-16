#!/usr/bin/env python3
# A2: run the pipeline's analyze_isoforms() (verbatim copy from DOOVER scripts/01_Description/isoform_NMD.py)
# on ALL genes of UniprotKB_isoforms_mapping_forNMD_noAmbiguous.bed (original run: DCR genes only).
# Only changes: input paths, gene list (chunk k of n), output path.
"""nmd_allgenes.py -- A2 step 1b: pipeline NMD calls (50-nt EJC rule) for the isoforms of ALL genes.

Purpose  analyze_isoforms() copied verbatim from DOOVER scripts/01_Description/isoform_NMD.py (the
         original run covered DCR genes only), run on chunk k of n of all genes.
Inputs   --bed  $DOOVER/output/01_Description/UniprotKB_isoforms_mapping_forNMD_noAmbiguous.bed
Outputs  --outdir/nmd_chunks/NMDtable_all.<k>.tsv ; concatenate the chunks into
         --outdir/NMDtable_allgenes.tsv (run_revision.sh does this).
         The verbatim function appends errors to ./output/01_Description/NMD_errors.err relative to
         the working directory; this wrapper chdirs into --outdir first.
Usage    python nmd_allgenes.py K N [--bed ...] [--outdir ...]     (bedtools on PATH, pybedtools)
Manuscript  Supplementary Figure S10 (NMD panel); docs/A2_results.md.
"""
import sys, os, argparse
import pandas as pd
from pybedtools import BedTool as bt
debug=False
_ap = argparse.ArgumentParser()
_ap.add_argument("k", type=int); _ap.add_argument("n", type=int)
_ap.add_argument("--bed", default=os.path.join(os.environ.get("DOOVER", "."), "output", "01_Description",
                                               "UniprotKB_isoforms_mapping_forNMD_noAmbiguous.bed"))
_ap.add_argument("--outdir", default=os.path.join(os.environ.get("WORK", "work"), "splicing_vs_AS"))
_a = _ap.parse_args()
k, n = _a.k, _a.n
_bed = os.path.abspath(_a.bed)
os.makedirs(os.path.join(_a.outdir, "nmd_chunks"), exist_ok=True)
os.chdir(_a.outdir)
os.makedirs("output/01_Description",exist_ok=True)   # error log path used inside the verbatim function
isoforms         = pd.read_csv(_bed, sep = "\t", header=None)
isoforms.columns = ["chr", "start", "end", "isoform", "strand", "gene"]
isoforms         = isoforms[["chr", "start", "end", "isoform", "gene", "strand"]]
allgenes         = sorted(pd.unique(isoforms['gene']))
genes            = allgenes[k::n]
# the verbatim function re-uses EJCindex from the previous gene when a gene model has a single merged exon
# (IndexError branch); a chunk must therefore not START with such a gene (would raise UnboundLocalError).
# Put the first gene with >=2 distinct exon starts at the front; downstream semantics are unchanged.
_ns = isoforms[isoforms['gene'].isin(set(genes))].groupby('gene')['start'].nunique()
_first = next(x for x in genes if _ns.get(x, 0) >= 2)
genes = [_first] + [x for x in genes if x != _first]
isoforms         = isoforms[isoforms['gene'].isin(set(genes))]
totalIsoformsIn  = str(len(pd.unique(isoforms['isoform'])))
totalGenesIn     = str(len(genes))
### functions #####################################################################################

def analyze_isoforms():
	"""
	For each gene, create a model by merging isoform
	For each isoform, predict chance to undergo NMD based on the gene model
	"""
	print("[analyze isoforms] started")
	print("[analyze isoforms] " + totalIsoformsIn + " isofoms from " + totalGenesIn + " genes in input dataset")
	g = 0 # gene counter
	total_g = len(genes)
	isoNMD = []
	for gene in genes:
		g += 1
		if g == 1:
			print("[analyze isoforms] " + str(g) + "/" + str(total_g), end = "")
		elif g > 1:
			print("\r[analyze isoforms] " + str(g) + "/" + str(total_g), end = "")
		# print("---")
		# print(gene)
		geneTable = isoforms[isoforms['gene'] == gene]
		#print(geneTable)
		geneIso   = pd.unique(geneTable['isoform'])
		geneBed   = bt.from_dataframe(geneTable.sort_values(by=['start','end']))
		
		# test if all exons in the gene are in the same orientation!!!
		if len(pd.unique(geneTable['strand'])) > 1:
			with open("output/01_Description/NMD_errors.err", "a") as errorfile:
				errorfile.write(str(gene) + " (conflicting exon directions in gene model)" + "\n")
			continue

		# merge the isoforms
		try:
			geneModel = geneBed.merge(s=True)
			#print(geneModel)
		except Exception as e:
			#print(" merge error for " + str(gene))
			#print(e)
			#print("\n")
			with open("output/01_Description/NMD_errors.err", "a") as errorfile:
				print("\n[analyze isoforms][warning]" + str(gene) + "(bedtools merge error [probably mapped to >1 chromosome]) gene skipped!\n")
				errorfile.write(str(gene) + "(bedtools merge error [probably mapped to >1 chromosome])\n")
			continue			

		for isoform in geneIso:
			isoTable  = geneTable[geneTable['isoform'] == isoform]
			isoStrand = pd.unique(isoTable['strand'])
			if debug:
				print("---")
				print(gene + " " + isoform + " " + isoStrand)
			isoExonsBed = bt.from_dataframe(isoTable.sort_values(by=['start','end']))
			geneModelXiso = geneModel.intersect(isoExonsBed, wao=True)
			# convert to list
			geneModelXisoL = [list(interval) for interval in geneModelXiso]
			# then to pandas df
			geneModelXisoDF = pd.DataFrame(geneModelXisoL)
			geneModelXisoDF.columns = ["ModelChrom", "ModelStart", "ModelEnd", "IsoChrom", "IsoStart", "IsoEnd", "isoform", "gene", "strand", "ol"]
			if debug:
				print(geneModelXisoDF)

			if isoStrand == "+":
				# find EJC on model
				try:
					modelStarts = sorted(pd.unique(geneModelXisoDF['ModelEnd'])) # get the list of unique exons start
					EJC = int(modelStarts[-2]) # one before last
					if debug:
						print("EJC = " + str(EJC))#)
					EJCindices = [i for i, x in enumerate(geneModelXisoDF['ModelEnd']) if int(x) == EJC]
					if debug:
						print(EJCindices)
					EJCindex = max(EJCindices)
					if debug:
						print("EJCindex = " + str(EJCindex))
				except IndexError:
					pass
				# find isoform end (potential PTC)
				PTC = max(pd.unique(geneModelXisoDF['IsoEnd']))
				if debug:
					print("PTC = " + str(PTC))
				# find PTC index in the intersect table
				PTCindices = [i for i, x in enumerate(geneModelXisoDF['IsoEnd']) if x == PTC]
				PTCindex = max(PTCindices)
				if PTCindex > EJCindex: # if the isoform ends after the last E-E junction
					# write no NMD
					NMDcall = "no_NMD"
					#print("No NMD; EJC-PTC distance = " + str(NMDdistance) + " bp")					
					isoNMD.append([gene, isoform, isoStrand, "last_exon", NMDcall])
				else:
					#print("PTCindex = " + str(PTCindex))
					# calculate the difference between the model pre-EJC exon and the PTC
					EJdiff = int(geneModelXisoDF['ModelEnd'][PTCindex]) - int(geneModelXisoDF['IsoEnd'][PTCindex])
					#print("EJdiff = " + str(EJdiff))
					# check if we need to add exons
					exonsBpToAdd = 0
					if abs(EJCindex - PTCindex) > 1:
						for exon in range(PTCindex + 1, EJCindex):
							try:
								exonLen = int(geneModelXisoDF['ModelEnd'][exon]) - int(geneModelXisoDF['ModelStart'][exon])# do sum of each exon in model
								exonsBpToAdd += exonLen
							except KeyError as ke: # it means we have 1 or 2 lines in the bed and thus there is nothing to add
								exonsBpToAdd = 0
								print("\n[analyze isoforms][warning]" + str(gene) + " " + isoform + " (no extra exon to add)" + "\n")
								with open("output/01_Description/NMD_errors.err", "a") as errorfile:
									errorfile.write(str(gene) + " " + isoform + " (no extra exon to add)" + "\n")
								pass
					#print("bp to add to EJdiff = " + str(exonsBpToAdd))
					NMDdistance = int(EJdiff) + int(exonsBpToAdd)
					if NMDdistance < 50:
						NMDcall = "no_NMD"
						#print("No NMD; EJC-PTC distance = " + str(NMDdistance) + " bp")					
						isoNMD.append([gene, isoform, isoStrand, NMDdistance, NMDcall])
					else:
						NMDcall = "NMD"
						#print("NMD; EJC-PTC distance = " + str(NMDdistance) + " bp")
						isoNMD.append([gene, isoform, isoStrand, NMDdistance, NMDcall])

			elif isoStrand == "-":
				try:
					modelStarts = sorted(pd.unique(geneModelXisoDF['ModelStart'])) # get the list of unique exons start (end column because in - strand)
					EJC = int(modelStarts[1]) + 1 # add + 1 to match coordinates system of +
					if debug:
						print("EJC = " + str(EJC))#)
					EJCindices = [i for i, x in enumerate(geneModelXisoDF['ModelStart']) if int(x) == (EJC-1)]
					EJCindex = min(EJCindices)
					if debug:
						print("EJCindex = " + str(EJCindex))				
				except IndexError:
					pass

				cleanIsoStarts = [num for num in geneModelXisoDF['IsoStart'] if int(num) >= 0] # IsoStarts are the end (-1) of the each exons
				##print(cleanIsoStarts)
				PTC = min(cleanIsoStarts)
				if debug:
					print("PTC = " + str(PTC))
				# find PTC index in the intersect table
				PTCindices = [i for i, x in enumerate(geneModelXisoDF['IsoStart']) if x == PTC] # t
				PTCindex = min(PTCindices)
				if PTCindex < EJCindex:
					# no NMD
					NMDcall = "no_NMD"
					#print("No NMD; EJC-PTC distance = " + str(NMDdistance) + " bp")					
					isoNMD.append([gene, isoform, isoStrand, "last_exon", NMDcall])
				else:
					#print("PTCindex = " + str(PTCindex))		
					# calculate the difference between the model pre-EJC exon and the PTC
					EJdiff = abs(int(geneModelXisoDF['ModelStart'][PTCindex]) - int(geneModelXisoDF['IsoStart'][PTCindex]))
					#print("EJdiff = " + str(EJdiff))
					# check if we need to add exons
					exonsBpToAdd = 0
					if abs(EJCindex - PTCindex) > 1:
						for exon in range(EJCindex, PTCindex - 1):
							try:
								exonLen = int(geneModelXisoDF['ModelEnd'][exon]) - int(geneModelXisoDF['ModelStart'][exon])# do sum of each exon in model
								exonsBpToAdd += exonLen
							except KeyError as ke: # it means we have 1 or 2 lines in the bed and thus there is nothing to add
								print("\n[analyze isoforms][warning]" + str(gene) + " " + isoform + " (no extra exon to add)" + "\n")
								exonsBpToAdd = 0
								with open("output/01_Description/NMD_errors.err", "a") as errorfile:
									errorfile.write(str(gene) + " " + isoform + " (no extra exon to add)" + "\n")
								pass
					#print("bp to add to EJdiff = " + str(exonsBpToAdd))
					NMDdistance = int(EJdiff) + int(exonsBpToAdd)
					if NMDdistance < 50:
						NMDcall = "no_NMD"
						#print("No NMD; EJC-PTC distance = " + str(NMDdistance) + " bp")					
						isoNMD.append([gene, isoform, isoStrand, NMDdistance, NMDcall])
					else:
						NMDcall = "NMD"
						#print("NMD; EJC-PTC distance = " + str(NMDdistance) + " bp")
						isoNMD.append([gene, isoform, isoStrand, NMDdistance, NMDcall])

	return(isoNMD)


rows=analyze_isoforms()
with open(f"nmd_chunks/NMDtable_all.{k}.tsv","w") as f:
    for row in rows:
        row=list(row); row[2]=str([str(x) for x in row[2]])  # serialise strand as in original NMDtable.tsv (pandas-version independent)
        f.write("\t".join(map(str,row))+"\n")
