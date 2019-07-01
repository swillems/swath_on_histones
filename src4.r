dataset = read.csv(
  '/media/NAS_proteomics/Sander/Projects/Laura/Relative_abundance/190129_HDAC_BCS_SWATHpaper/190130_HDAC_RA_SW.csv',
  header=TRUE,
  stringsAsFactors=FALSE,
  sep=";"
)

attach(dataset)

Time = as.factor(Time)

outfile_name = '/media/NAS_proteomics/Sander/Projects/Laura/Relative_abundance/190129_HDAC_BCS_SWATHpaper/190130_HDAC_RA_SW_parsed.csv'

write.table("", outfile_name, sep=",", col.names=FALSE, row.names=FALSE)

# peptide = GKGGKGLGKGGAKR_.5..EG.Succinyl..K.
for (peptide_name in colnames(dataset[2:dim(dataset)[2]])){
  # plot(as.numeric(Time), peptide)
  peptide = unlist(dataset[peptide_name])
  fit = aov(peptide ~ Time)
  z = summary(fit)
  # write.table(z[[i]], outfile_name, sep=",", append=TRUE)
  z
  p = z[[1]]$`Pr(>F)`[1]
  # peptide_pvals[i] = p
  x = pairwise.t.test(peptide, Time, p.adj = "none")
  x
  y = x$p.value
  write.table(peptide_name, outfile_name, sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  write.table(t(c("Pr(>F)", p)), outfile_name, sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  # write.table(p, outfile_name, sep=",", append=TRUE, col.names = FALSE, row.names=FALSE)
  write.table(t(c("",colnames(y))), outfile_name, sep=",", append=TRUE, col.names=FALSE, row.names=FALSE)
  write.table(y, outfile_name, sep=",", append=TRUE, col.names=FALSE)
  write.table("", outfile_name, sep=",", append=TRUE, col.names = FALSE, row.names=FALSE)
  # peptide_pvals[i] = y[4, 1]
}

