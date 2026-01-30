# DCR-Project

Raw data and code supporting the Dual Coding Region publication (DOI TBD).

## Data

The `Data` directory contains the essential raw UniprotKB/SwissProt mapping of human protein isoforms performed wit MIRAGE2 on the human reference genome hg38.p14 as well as supplementary pre-formatted data and metadata required to reproduce the analyses presented in the paper.

| Filename | Description |
| -------- | ----------- |
| `Species-MSAs.tar.gz` | This archive contains the raw MIRAGE2 alignments of human UniprotKB/SwissProt isoforms to the human reference genome hg38.p14 which content is required by some scripts in the `Code` directory |
| `Extracted-DCEs_06202025.tar.gz` | This archive contains the raw extracted DCRs from the MIRAGE2 alignments. It includes the main tables and files consumed by the analysis pipeline found in the `Code` directory of this repository | 
| `demultiplexed.tallyman.counts.unique.wHeader.tsv` | GTEx tissue-specific DCR exon-exon junction counts from Tallyman (raw counts) |
| `GTExMetadataUpdated.tsv` | Metadata for the GTEx dataset |
| `Mice.tar.gz` | This archive contains the Mice expression data, including tissue-specific DCR exon-exon junction counts with Tallyman and metadata |
| `CIFs.tar.gz` | Alfafold3 structure predictions for DCRs (CIF format) |
| `connected_DCR_representatives.tsv` | Selection of isoforms pairs representing canonical and alternative DCR peptides to be used for structure prediction with Alfafold3 |
| `Table_S1.tsv` | Retained DCR with isoform metadata, including canonical/alternative frame adjudication |

## Code

The code is organized in 4 main sections with the main pipeline scripts located in the top directory level and dependent scripts in their respective sub-directories.

| Filename | Description |
| -------- | ----------- |
| **`01_Description.sh`** | Process the raw data, filter DCRs and generate main statistics. Produce outputs used by all subsequent analysis steps and panels for **Figure 3** |
| **`02_Conservation.sh`** | Perform data analyses related to conservation analyses between human and mice and PhyloP conservation. Produce panels for **Figure 2** |
| **`03_GTEx_part1.sh`** | Format the data for RNA-seq quantification with Tallyman |
| **`03_GTEx_part2.sh`** | Handle post quantification data, produce figure S5 and launch the differential frame usage analysis. |
| `03_GTEx/GTEx_Plot_DotMatrix.R` | Figure 4 |
| `03_GTEx/Mice.R` | Figure 5 |
| `04_Peptide_function/01_extract_isoforms.py` | Extract the DCR peptides in fasta format |
| `04_Peptide_function/02_fasta_to_afserver_batches.py` | Creates Alfafold3 server batch files |
| `04_Peptide_function/03_atom2res_plddt.py` | Extract per atom plddts to AA for Alfafold3 CIF files |
| `04_Peptide_function/04_extract_DCR_plddts.sh` | Format plddts |
| `04_Peptide_function/05_batch_dssp.py` | Run DSSP in batch over the CIF files to extract structures |
| `04_Peptide_function/06_parse_dssp.sh` | Format DSSP structures |
|`04_Peptide_function/07_plot_AF3_DCR_plddts.R` | Figure 6 |

