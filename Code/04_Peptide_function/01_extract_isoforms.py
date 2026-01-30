#!/usr/bin/env python3
import pandas as pd
import re
from Bio import SeqIO
import os

input_table = "Data/connected_DCR_representatives.tsv"
canonical_fasta = "Data/uniprot_sprot.fasta" # to obtain from UniprotKB/Swiss-Prot
alt_fasta = "Data/uniprot_sprot_varsplic.fasta" # to obtain from UniprotKB/Swiss-Prot

DCR_list = pd.read_csv(input_table, sep="\t", header=0, dtype=str)

output_dir = "output/04_Peptide_function/"
os.makedirs(output_dir, exist_ok=True)

def extract_canonical(DCR_list):
    """
    Extract canonical sequences from the UniProt FASTA files based on DCR IDs.
    """

    for DCR_ID in DCR_list['DCR']:
        #print(f"Processing DCR ID: {DCR_ID}")
        can_iso = DCR_list[DCR_list['DCR'] == DCR_ID]['canonical'].values[0]
        prot_id = DCR_list[DCR_list['DCR'] == DCR_ID]['uniprot'].values[0]
        gene = DCR_list[DCR_list['DCR'] == DCR_ID]['gene_regex'].values[0]
        frame = DCR_list[DCR_list['DCR'] == DCR_ID]['can_frame'].values[0]
        print(f"Canonical isoform for {DCR_ID}: {can_iso}; {prot_id}; {gene}")
        if can_iso == "1":
            #print("searching in canonical fasta")
            with open(canonical_fasta, "r") as fasta:
                for line in fasta:
                    if line.startswith(">"):
                        if prot_id in line:
                            header_line = line.strip()
                            header = header_line.split(' ')[0].replace('>', '')
                            #print(f"Found header line: {header_line}")
                            #print(f"Processed header: {header}")
                            for record in SeqIO.parse(canonical_fasta, "fasta"):
                                if record.id == header:
                                    #print(f"Found sequence for {header}")
                                    #print(record.seq)
                                    new_header = f"{prot_id}-{can_iso}|{gene}|DCR={DCR_ID}|frame={frame}|canonical"
                                    record.id = new_header
                                    record.description = ""
                                    output_filename = f"{output_dir}{DCR_ID}.{gene}.fasta"
                                    with open(output_filename, "w") as out_fasta:
                                        SeqIO.write(record, out_fasta, "fasta")
                                    #print(record.seq)
                                    break
        else:
            with open(alt_fasta, "r") as fasta:
                for line in fasta:
                    if line.startswith(">"):
                        iso_id = prot_id + "-" + can_iso
                        if iso_id in line:
                            header_line = line.strip()
                            header = header_line.split(' ')[0].replace('>', '')
                            #print(f"Found header line: {header_line}")
                            #print(f"Processed header: {header}")
                            for record in SeqIO.parse(alt_fasta, "fasta"):
                                if record.id == header:
                                    #print(f"Found sequence for {header}")
                                    #print(record.seq)
                                    new_header = f"{prot_id}-{can_iso}|{gene}|DCR={DCR_ID}|canonical"
                                    record.id = new_header
                                    record.description = ""
                                    output_filename = f"{output_dir}{DCR_ID}.{gene}.fasta"
                                    with open(output_filename, "w") as out_fasta:
                                        SeqIO.write(record, out_fasta, "fasta")
                                    #print(record.seq)
                                    break

def extract_alternative(DCR_list):
    """
    Extract alternative sequences from the UniProt FASTA files based on DCR IDs.
    """

    for DCR_ID in DCR_list['DCR']:
        can_iso = DCR_list[DCR_list['DCR'] == DCR_ID]['alt'].values[0]
        prot_id = DCR_list[DCR_list['DCR'] == DCR_ID]['uniprot'].values[0]
        gene = DCR_list[DCR_list['DCR'] == DCR_ID]['gene_regex'].values[0]
        frame = DCR_list[DCR_list['DCR'] == DCR_ID]['alt_frame'].values[0]
        print(f"Alternative isoform for {DCR_ID}: {can_iso}; {prot_id}; {gene}")
        with open(alt_fasta, "r") as fasta:
            for line in fasta:
                if line.startswith(">"):
                    iso_id = prot_id + "-" + can_iso
                    if iso_id in line:
                        header_line = line.strip()
                        header = header_line.split(' ')[0].replace('>', '')
                        for record in SeqIO.parse(alt_fasta, "fasta"):
                            if record.id == header:
                                if can_iso == "1":
                                    new_header = f"{prot_id}-{can_iso}.1|{gene}|DCR={DCR_ID}|frame={frame}|alt"
                                else:
                                    new_header = f"{prot_id}-{can_iso}|{gene}|DCR={DCR_ID}|frame={frame}|alt"
                                record.id = new_header
                                record.description = ""
                                output_filename = f"{output_dir}{DCR_ID}.{gene}.fasta"
                                with open(output_filename, "a") as out_fasta:
                                    SeqIO.write(record, out_fasta, "fasta")
                                break

extract_canonical(DCR_list)
extract_alternative(DCR_list)