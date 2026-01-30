# DCR plddts 
library(ggplot2)
library(dplyr)
library(reshape2)
library(ggsankey)

plddts<-read.table("extract_DCR_plddts.tsv")
names(plddts)<-c("DCR_ID", "gene", "isoform", "PepPos", "AA", "plddt", "form")

# control canonical/alternative status since update
updated_frame_assignment<-read.table("Table_S1.tsv", h = T)
pep_iso_slice <- updated_frame_assignment[updated_frame_assignment$DCR %in% unique(plddts$DCR_ID),] 
pep_iso_slice %>% separate(col = isoform, into = c("species", "gene", "iso_number"), sep = "\\|") -> pep_iso_slice_processed
pep_iso_slice_processed2<-pep_iso_slice_processed[,c("DCR", "iso_number", "new_frame_status")]
pep_iso_slice_processed2$DCR_iso <- paste(pep_iso_slice_processed2$DCR, pep_iso_slice_processed2$iso_number, sep = "_")

plddts_summary <- plddts %>% group_by(DCR_ID, isoform) %>% summarise(form = unique(form))
plddts_summary$DCR_iso <- paste(plddts_summary$DCR_ID, plddts_summary$isoform, sep = "_")

pep_iso_slice_processed2 <- pep_iso_slice_processed2[pep_iso_slice_processed2$DCR_iso %in% unique(plddts_summary$DCR_iso),]

comparison <- merge(plddts_summary, pep_iso_slice_processed2, by = "DCR_iso")


# structure level
dssp<-read.table("DCR_peps_structure_annotation.tsv")
names(dssp)<-c("DCR", "isoform", "frame", "pos", "AA", "structure", "form")

perIsoDCR<-unique(updated_frame_assignment[,c("DCR", "frame", "new_frame_status")])
perIsoDCR$key<-paste(perIsoDCR$DCR, perIsoDCR$frame, sep = "_")
perIsoDCR$frame_status <- perIsoDCR$new_frame_status
iso_form<-unique(perIsoDCR[c("key", "frame_status")])
dssp$key<-paste(dssp$DCR, dssp$frame, sep = "_")
dssp<-merge(dssp[,c(1:6,8)], iso_form, by = "key")

plddts$AA_key <- paste(plddts$DCR_ID, plddts$isoform, plddts$PepPos, plddts$AA, sep = "_")

dssp[dssp$isoform == "1.1",]$isoform <- "1-1"
dssp$AA_key <- paste(dssp$DCR, dssp$isoform, dssp$pos, dssp$AA, sep = "_")
struct_plddt<-merge(plddts[,c("AA_key", "plddt")], dssp, by = "AA_key")
struct_plddt$confident_structure <- ifelse(struct_plddt$plddt >= 70, struct_plddt$structure, "low_plddt")

struct_plddt2<-struct_plddt %>% group_by(DCR, frame_status, confident_structure) %>% summarize(struc_count = length(confident_structure))
struct_plddt2$frame_status<-factor(struct_plddt2$frame_status, levels = c("canonical", "alternative"))

struct_plddt2$confident_structure <- factor(struct_plddt2$confident_structure, 
                                            levels = c("low_plddt", ".", "H", "B", "E", "G", "I", "P", "T", "S"),
                                            labels = c("AA plddt < 70", "disordered", "α-helix", "residue in isolated β-bridge", 
                                                       "extended strand, participates in β ladder", "3[10]-helix", "π-helix", "poly-proline II helix",
                                                       "hydrogen-bonded turn", "bend"))
ggplot(struct_plddt2)+
  geom_bar(aes(x = frame_status, y = struc_count, fill = confident_structure), position = "fill", stat = "identity")+
  ylab("proportion among isoforms")+
  xlab("isoform DCR category")+
  scale_fill_manual(name = "structure at AA",
                    values = c("AA plddt < 70" = "grey80",
                               "disordered" = "grey30",
                               "α-helix" = "#4e79a7",
                               "residue in isolated β-bridge" = "#edc948", 
                               "extended strand, participates in β ladder" = "#f28e2b", 
                               "3[10]-helix" = "#59a14f",
                               "π-helix" = "#76b7b2", 
                               "poly-proline II helix" = "#b07aa1", 
                               "hydrogen-bonded turn" = "#ff9da7",
                               "bend" = "#e15759"),
                    labels = c("AA plddt < 70",
                               "disordered",
                               "α-helix",
                               "residue in isolated β-bridge",
                               "extended strand, participates in β ladder",
                               expression(3[10]-helix),
                               "π-helix",
                               "poly-proline II helix",
                               "hydrogen-bonded turn",
                               "bend")
                    )+
  theme_classic()+
  theme(legend.position = "bottom", legend.direction = "vertical", axis.text = element_text(size = 14), axis.title = element_text(size = 18), legend.text = element_text(size = 12), legend.title = element_text(size = 16))+
  guides(fill = guide_legend(ncol = 2))

### stats

length(unique(struct_plddt2$DCR))
