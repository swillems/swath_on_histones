# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("limma")

#input = "HDAC"
#input = "HDAC_norm"
#input = "BCC_453"
input = "BCC_468"

if(input=="HDAC"){
  infile_name = '/home/sander/Documents/LabFBT/Colleague_work/180306_LdC_SWATH/Relative_abundance/190129_HDAC_BCS_SWATHpaper/190130_HDAC_RA_SW.csv'
  outfile_name = '/home/sander/Documents/LabFBT/Colleague_work/180306_LdC_SWATH/Relative_abundance/190129_HDAC_BCS_SWATHpaper/190130_HDAC_RA_SW_limma.csv'
  all_contrasts = list(
    "Time1-Time2",
    "Time1-Time3",
    "Time1-Time4",
    "Time1-Time5",
    "Time1-Time6",
    "Time1-Time7",
    "Time1-Time8",
    "Time1-Time9",
    "Time2-Time3",
    "Time2-Time4",
    "Time2-Time5",
    "Time2-Time6",
    "Time2-Time7",
    "Time2-Time8",
    "Time2-Time9",
    "Time3-Time4",
    "Time3-Time5",
    "Time3-Time6",
    "Time3-Time7",
    "Time3-Time8",
    "Time3-Time9",
    "Time4-Time5",
    "Time4-Time6",
    "Time4-Time7",
    "Time4-Time8",
    "Time4-Time9",
    "Time5-Time6",
    "Time5-Time7",
    "Time5-Time8",
    "Time5-Time9",
    "Time6-Time7",
    "Time6-Time8",
    "Time6-Time9",
    "Time7-Time8",
    "Time7-Time9",
    "Time8-Time9"
  )
}
if(input=="HDAC_norm"){
  infile_name = '/home/sander/Documents/LabFBT/Colleague_work/180306_LdC_SWATH/Relative_abundance/190129_HDAC_BCS_SWATHpaper/181212_HDAC_norm_parse_edit_RA_forStats_SW.csv'
  outfile_name = '/home/sander/Documents/LabFBT/Colleague_work/180306_LdC_SWATH/Relative_abundance/190129_HDAC_BCS_SWATHpaper/181212_HDAC_norm_parse_edit_RA_forStats_SW_limma.csv'
  all_contrasts = list(
    "Time1-Time2",
    "Time1-Time3",
    "Time1-Time4",
    "Time1-Time5",
    "Time1-Time6",
    "Time1-Time7",
    "Time1-Time8",
    "Time1-Time9",
    "Time2-Time3",
    "Time2-Time4",
    "Time2-Time5",
    "Time2-Time6",
    "Time2-Time7",
    "Time2-Time8",
    "Time2-Time9",
    "Time3-Time4",
    "Time3-Time5",
    "Time3-Time6",
    "Time3-Time7",
    "Time3-Time8",
    "Time3-Time9",
    "Time4-Time5",
    "Time4-Time6",
    "Time4-Time7",
    "Time4-Time8",
    "Time4-Time9",
    "Time5-Time6",
    "Time5-Time7",
    "Time5-Time8",
    "Time5-Time9",
    "Time6-Time7",
    "Time6-Time8",
    "Time6-Time9",
    "Time7-Time8",
    "Time7-Time9",
    "Time8-Time9"
  )
}
if(input=="BCC_453"){
  infile_name = '/home/sander/Documents/LabFBT/Colleague_work/180306_LdC_SWATH/Relative_abundance/190129_HDAC_BCS_SWATHpaper/181211_BCC_Norm_parse_edit_RA_forStats_453_SW.csv'
  outfile_name = '/home/sander/Documents/LabFBT/Colleague_work/180306_LdC_SWATH/Relative_abundance/190129_HDAC_BCS_SWATHpaper/181211_BCC_Norm_parse_edit_RA_forStats_453_SW_limma.csv'
  all_contrasts = list(
    "Time1-Time2",
    "Time1-Time3",
    "Time2-Time3"
  )
}
if(input=="BCC_468"){
  infile_name = '/home/sander/Documents/LabFBT/Colleague_work/180306_LdC_SWATH/Relative_abundance/190129_HDAC_BCS_SWATHpaper/181211_BCC_Norm_parse_edit_RA_forStats_468_SW.csv'
  outfile_name = '/home/sander/Documents/LabFBT/Colleague_work/180306_LdC_SWATH/Relative_abundance/190129_HDAC_BCS_SWATHpaper/181211_BCC_Norm_parse_edit_RA_forStats_468_SW_limma.csv'
  all_contrasts = list(
    "Time1-Time2",
    "Time1-Time3",
    "Time2-Time3"
  )
}


dataset = read.csv(
  infile_name,
  header=TRUE,
  stringsAsFactors=FALSE,
  sep=";"
)

write.table(t(colnames(dataset[2:dim(dataset)[2]])), outfile_name, sep=",", col.names=FALSE)

library(limma)
# transpose of dataset 
dataset_t <- t(dataset)
colnames(dataset_t) <- dataset_t[1, ]
dataset_t <- dataset_t[-1, ]

# Create design matrix
target <- data.frame(Time = colnames(dataset_t))
target$Time <- as.factor(target$Time)

design <- model.matrix(~ 0+ Time, data = target)

# Fit full model 
fit <- lmFit(dataset_t, design = design)


# Look at contrasts
for(current_contrast in all_contrasts){
  write.table(current_contrast, outfile_name, sep=",", append=TRUE, col.names = FALSE, row.names=FALSE)
  contrast <- makeContrasts(current_contrast, levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  z = topTable(fit2, adjust="BH", number = Inf, sort.by = "none")
  write.table(t(z), outfile_name, sep=",", append=TRUE, col.names=FALSE)
}

