# the below piece of code was created to make complete summary pictures of samples
# they include the error plots of the samples with and without penalty
# the code is easily adjusted to test and plot other variables (different ploidy, different error models, etc)
# the variable models refers to a dataframe of a completed models.txt file

for (a in 1:length(pd$name)) {
  error0 <- singlemodel(matias96,QDNAseqobjectsample = a,penalty = 0)
  error0.5 <- singlemodel(matias96,QDNAseqobjectsample = a,penalty = 0.5)
  min0 <- error0$minima[which(error0$rerror==min(error0$rerror))]
  last0 <- tail(error0$minima,1)
  min0.5 <- error0.5$minima[which(error0.5$rerror==min(error0.5$rerror))]
  last0.5 <- tail(error0.5$minima,1)
  list <- unique(c(min0,min0.5,last0,last0.5))
  plots <- vector(mode = 'list', length = 3+length(list))
  plots[[1]] <- error0$errorplot + ggtitle("errorlist - penalty = 0")
  plots[[2]] <- error0.5$errorplot + ggtitle("errorlist - penalty = 0.5")
  st <- models$standard[which(models$sample==pd$name[a])]
  pl <- models$ploidy[which(models$sample==pd$name[a])]
  cell <- models$likely_fit[which(models$sample==pd$name[a])]
  plots[[3]] <- singleplot(matias96,QDNAseqobjectsample = a,standard = st, ploidy = pl, cellularity = cell, title = paste0(pd$name[a], " - chosen model"))
  for (m in 1:length(list)) {
    plots[[3+m]] <- singleplot(matias96,QDNAseqobjectsample = a, cellularity = list[m], title = pd$name[a])
  }
  dir.create("all_plots")
  if(length(plots)==7) {
    png(file.path("all_plots",paste0(pd$name[a],".png")),width = 1920, height = 960)
    print(multiplot(plotlist = plots, layout = matrix(c(1,2,"",3,4,5,6,7), nrow=2, byrow = TRUE)))
    dev.off()
  } else {
    png(file.path("all_plots",paste0(pd$name[a],".png")),width = 1440, height = 960)
    print(multiplot(plotlist = plots, cols=3))
    dev.off()
  }
}
