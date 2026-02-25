library(UpSetR)

# Define your sets
DCRdata <- read.table("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/output/01_Description/filtered_homo_sapiens_DCR.tsv", h = T)
Alex <- toupper(unique(DCRdata$Gene))
# Other study sets
Chung2007 <- readLines("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/data/paper_formatted_gene_lists/Chung2007/Chung2007.uniq.txt")
Kovaks2010 <- readLines("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/data/paper_formatted_gene_lists/Kovaks2010/kovaks2010_formatted.txt")
Liang2006 <- readLines("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/data/paper_formatted_gene_lists/Liang2006/Liang2006_INCOMPLETE_.uniq.txt")
Lin2011 <- readLines("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/data/paper_formatted_gene_lists/Lin2011/Lin2011.uniq.txt")
Michel2012 <- readLines("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/data/paper_formatted_gene_lists/Michel2012/Michel2012.uniq.txt")
Ribrioux2008 <- readLines("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/data/paper_formatted_gene_lists/Ribrioux2008/Ribrioux2008_INCOMPLETE_.uniq.txt")
VDP2013mRNA <- readLines("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/data/paper_formatted_gene_lists/Vanderperre2013/vanderperre2013.mRNA.allpredictions.uniq.txt")
VDP2013prot <- readLines("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/data/paper_formatted_gene_lists/Vanderperre2013/vanderperre2013.PROTEINS.uniq.alldatasets.txt")
Xu2010 <- readLines("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2026/02_04_DCR-merge-overlap/DOOVER/data/paper_formatted_gene_lists/Xu2010/Xu2010.uniq.txt")

# Convert to UpSet format
sets_list <- list(ThisStudy = Alex, Chung2007 = Chung2007, Kovacs2010 = Kovaks2010, Liang2006 = Liang2006, Lin2011 = Lin2011,
                  Michel2012 = Michel2012, Ribrioux2008 = Ribrioux2008, VDP2013mRNA = VDP2013mRNA, VDP2013prot = VDP2013prot, 
                  Xu2010 = Xu2010)
sets_matrix <- fromList(sets_list)

# Find all sets with an interesection with "ThisStudy"
rSums = rowSums(sets_matrix)

# Find observed intersections
observed_intersections <- list()
for (i in 1:nrow(sets_matrix)) {
  active_sets <- names(sets_matrix)[sets_matrix[i, ] == 1]
  observed_intersections <- append(observed_intersections, list(active_sets))
}

# Remove duplicates and keep only sets containing "ThisStudy"
valid_intersections <- unique(observed_intersections)
valid_intersections <- valid_intersections[sapply(valid_intersections, function(x) "ThisStudy" %in% x)]

# Create UpSetR queries dynamically
qs <- lapply(valid_intersections, function(intersection) {
  list(query = intersects, params = intersection, color = "orange", active = TRUE)
})

# Generate UpSet plot (for supp mat probably)
upset(sets_matrix, nintersects = NA,
      sets = names(sets_list),
      order.by = "freq", 
      mainbar.y.label = "DCR Genes \n(intersection size)",
      queries = qs, number.angles = 0,
      sets.x.label = "Gene Set Size",
      sets.bar.color = c("blue4", "blue4", "#aa72c4", "#aa72c4", "blue4","#aa72c4", "blue4", "blue4", "blue4", "#e3a9e3")
)

length(Alex)

allStudies<-unique(c(Chung2007, Kovaks2010, Liang2006, Lin2011, Michel2012, Ribrioux2008, VDP2013mRNA, VDP2013prot, Xu2010))
length(Alex[!Alex %in% allStudies])
length(Alex[Alex %in% VDP2013mRNA])
