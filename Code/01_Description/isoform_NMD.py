#!/usr/local/bin/python
import pandas as pd
from pybedtools import BedTool as bt
import os
import argparse

### parser ########################################################################################

parser = argparse.ArgumentParser(description="Predict NMD from UniprotKB isoforms and Mirage DCRs")
parser.add_argument("-n", "--nmdtable", action="store_true", help="starts from NMD table (step2)")
parser.add_argument("-d", "--debug", action="store_true", help="becomes talkative")
# parse args
args = parser.parse_args()
fromNMDtable = args.nmdtable
debug = args.debug

### import data and initialize variables ##########################################################

DCR              = pd.read_csv("output/01_Description/filtered_DCR_overlapping_regions_with_gene_and_dceID.bed", sep = "\t", header=None)
DCR.columns      = ["chr", "start", "end", "DCRindex", "strand", "gene", "frame"]
isoforms         = pd.read_csv("output/01_Description/UniprotKB_isoforms_mapping_forNMD_noAmbiguous.bed", sep = "\t", header=None)
isoforms.columns = ["chr", "start", "end", "isoform", "strand", "gene"]
isoforms         = isoforms[["chr", "start", "end", "isoform", "gene", "strand"]]
totalIsoformsIn  = str(len(pd.unique(isoforms['isoform'])))
genes            = pd.unique(DCR['gene'])
#genes            = ["atp2b3"]
totalGenesIn     = str(len(genes))
DCRiso           = pd.read_csv("output/01_Description/isoforms_to_DCR_metadata.tsv", sep = "\t", header=None)
DCRiso.columns   = ['isoform', 'frame', 'termini', 'DCR']
# clean error file
if fromNMDtable == False:
	try:
		os.remove("NMD_errors.err")
	except FileNotFoundError:
		pass

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

def analyze_DCR(NMDtable, DCRiso):
	"""
	Assess for each DCR case if all isoform of one frame goes NMD
	This means that the DCR leads to NMD
	However, if one frame is all NMD but the termini is also found in anoter frame
	that is not 100% NMD, the NMD is NOT caused by DCR
	"""
	DCRisoTable = NMDtable.merge(DCRiso, on='isoform', how='inner')
	print("[analyze DCR] after merging, " + str(len(pd.unique(DCRisoTable['isoform']))) + "/" + totalIsoformsIn + " isoforms are left (" + str(len(pd.unique(DCRisoTable['gene']))) + "/" + totalGenesIn + " genes; " + str(len(pd.unique(DCRisoTable['DCR']))) + " DCRs)" )
	DCRlist     = pd.unique(DCRisoTable['DCR'])
	DCRtot      = len(DCRlist)
	results = [] # create the list that will store the results
	DCRcounter = 0
	for DCR in DCRlist:
		DCRcounter += 1 # increment
		if DCRcounter == 1:
			print("[analyze DCR] " + str(DCRcounter) + "/" + str(DCRtot), end = "")
		elif DCRcounter > 1:
			print("\r[analyze DCR] " + str(DCRcounter) + "/" + str(DCRtot), end = "")
		perDCRtable = DCRisoTable[DCRisoTable['DCR'] == DCR]
		frames      = pd.unique(perDCRtable['frame'])
		frNMD       = [] # will store % NMD per frame
		if debug:
			print("---")
			print(DCR)
			print(perDCRtable)
		# Count occurrences of "NMD" and "no_NMD" per frame  
		counts = perDCRtable.groupby(["frame", "pred"]).size().unstack(fill_value=0)  
		# Compute ratio NMD / (NMD + no_NMD)  
		counts["ratio"] = counts.get("NMD", 0) / (counts.get("NMD", 0) + counts.get("no_NMD", 0))  
		if debug:
			print(counts)
		# Get frames where ratio == 1
		frames_with_ratio_1 = counts[counts["ratio"] == 1].index.tolist()

		# Check conditions  
		has_ratio_1 = len(frames_with_ratio_1) > 0
		has_ratio_less_than_1 = (counts["ratio"] < 1).any()
		nonNMDframesTermini = perDCRtable[perDCRtable["pred"] == "no_NMD"]['termini'].unique()
		NMDframesTermini = perDCRtable[perDCRtable["frame"].isin(frames_with_ratio_1)]['termini'].unique()
		if debug:
			print("non NMD frames termini: " + str(nonNMDframesTermini))
			print("NMD frames termini: " + str(NMDframesTermini))
		#print(counts.columns)
		# Final evaluation
		#alltermini = NMDframesTermini.tolist().append(nonNMDframesTermini.tolist())
		NMDframesTerminiL = NMDframesTermini.tolist()
		#print(NMDframesTerminiL)
		nonNMDframesTerminiL = nonNMDframesTermini.tolist()
		alltermini = NMDframesTerminiL + nonNMDframesTerminiL
		#print(alltermini)
		if has_ratio_1 and has_ratio_less_than_1:
			if NMDframesTermini.all() != 'C-TERM': #and set(NMDframesTermini).intersection(set(nonNMDframesTermini)):
				results.append([DCR, 'no_NMD', 'not_C-terminal', frames_with_ratio_1])
				if debug:
					print("DCR do not cause NMD: other frames with no NMD have this termini")
			else:
				results.append([DCR, 'NMD', 'in_NMD_zone', frames_with_ratio_1])
				if debug:
				    print("DCR causes NMD: Frames with ratio == 1:", frames_with_ratio_1)
		else:  
			if debug:
				print("DCR do not cause NMD")
			if nonNMDframesTermini.all() == 'C-TERM': # and NMDframesTermini.all() == 'C-TERM':
				results.append([DCR, 'no_NMD', 'last_exon', frames_with_ratio_1])
			#elif len(NMDframesTermini) > 1 and any(NMDframesTermini) != 'C-TERM':
			elif all(end != 'C-TERM' for end in alltermini):
				results.append([DCR, 'no_NMD', 'not_C-terminal', frames_with_ratio_1])
			else:
				results.append([DCR, 'no_NMD', 'coding_long_enough', frames_with_ratio_1])
		# if debug:
		# 	print(results)	
	outputDF = pd.DataFrame(results)
	print("[analyze DCR] At then end, we have " + str(len(pd.unique(outputDF[0]))) + "/" + str(len(pd.unique(DCRisoTable['DCR']))) + " DCRs with results")
	return(outputDF)




### workflow ######################################################################################

# we first check if we need to run the first function
if fromNMDtable:
	# Read NMDtable from previous run Pandas DataFrame
	try:
		NMDtable = pd.read_csv("output/01_Description/NMDtable.tsv", sep="\t", header=None)
		print("[info] NMD table loaded")
	except FileNotFoundError:
		print("[error] File NMDtable.tsv missing")
else:
	NMDrawTable = analyze_isoforms()
	# Write to TSV (checkpoint)
	with open("output/01_Description/NMDtable.tsv", "w", newline="") as f:
	    for row in NMDrawTable:
	        f.write("\t".join(map(str, row)) + "\n")
	    print("\n[info] NMD table written to NMDtable.tsv")
	NMDtable = pd.DataFrame(NMDrawTable)

# then we look for DCR-induced NMD
NMDtable.columns = ['gene', 'isoform', 'strand', 'dist', 'pred']
outputDF = analyze_DCR(NMDtable, DCRiso)
# check if there was a previous output
try:
	os.remove("output/01_Description/DCR_NMD_analysis.out")
except FileNotFoundError:
	pass
# write new output file
# with open("DCR_NMD_analysis.out", "w", newline="") as o:
# 	for row in outputDF:
# 		o.write("\t".join(map(str, row)) + "\n")
outputDF.to_csv("output/01_Description/DCR_NMD_analysis.out", index = False, sep = "\t")
print("[analyze DCR] Completed")
if debug:
	print(outputDF)



