# DCR-Project

Raw data and code supporting the Dual Coding Region publication (BioRxiv: https://doi.org/10.64898/2026.02.24.707805).

## Data

The `Data/` directory contains the essential raw UniprotKB/Swiss-Prot mapping of human protein isoforms performed wit [`MIRAGE2`](https://github.com/TravisWheelerLab/Mirage) on the human reference genome [`hg38.p14`](https://hgdownload.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz) as well as supplementary pre-formatted data and metadata required to reproduce the analyses presented in the paper.

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

The `Code/` directory is organized in 4 main sections. The main pipeline scripts are located in the top level and dependent scripts in their respective sub-directories.

| Filename | Description |
| -------- | ----------- |
| **`01_Description.sh`** | Process the raw data, filter DCRs and generate main statistics. Produce outputs used by all subsequent analysis steps and panels for **Figure 3** |
| **`02_Conservation.sh`** | Perform data analyses related to conservation analyses between human and mice and PhyloP conservation. Produce panels for **Figure 2** |
| **`03_GTEx_part1.sh`** | Format the data for RNA-seq quantification with Tallyman |
| **`03_GTEx_part2.sh`** | Handle post quantification data, produce figure S5 and launch the differential frame usage analysis. |
| `03_GTEx/GTEx_Plot_DotMatrix.R` | **Figure 4** |
| `03_GTEx/Mice.R` | **Figure 5** |
| `04_Peptide_function/01_extract_isoforms.py` | Extract the DCR peptides in fasta format |
| `04_Peptide_function/02_fasta_to_afserver_batches.py` | Creates Alfafold3 server batch files |
| `04_Peptide_function/03_atom2res_plddt.py` | Extract per atom plddts to AA for Alfafold3 CIF files |
| `04_Peptide_function/04_extract_DCR_plddts.sh` | Format plddts |
| `04_Peptide_function/05_batch_dssp.py` | Run DSSP in batch over the CIF files to extract structures |
| `04_Peptide_function/06_parse_dssp.sh` | Format DSSP structures |
|`04_Peptide_function/07_plot_AF3_DCR_plddts.R` | **Figure 6** |

### Revision analyses (2026)

`Code/05_Revision/` holds the code added for the GBE revision: the alternate-frame selection test (**Figure 2B(iii)**), the comparison of DCRs with ordinary alternative splicing (**Figure S10**), the splice-site, expression-class and tissue-specificity controls (**Figure S9**), the GABBR1 and GNAS gene models (**Figures S6, S7**), the GC-content controls (**Figure S8**), and an archive of the earlier revision code, kept as executed. A single driver, `Code/05_Revision/run_revision.sh`, runs these analyses on the outputs of the pipeline above. See [`Code/05_Revision/README.md`](Code/05_Revision/README.md).


## Questions/Requests

Questions and requests can be send through the [issue tab](https://github.com/TravisWheelerLab/DCR-project/issues), or via email using the address for correspondence listed on the publication (https://doi.org/10.64898/2026.02.24.707805)