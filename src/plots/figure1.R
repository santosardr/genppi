library(ggplot2)
data=read.csv("pp.csv",header=TRUE,stringsAsFactors=FALSE, colClasses =c("character", "integer", "integer", "integer"));
pdf("Figure1-BoxPlot-PP-by-Genomes.pdf", width=11.69, height=8.27, pointsize=8, );
p<-ggplot(data,aes(x=Genome, y=Genomes))+geom_boxplot(varwidth = TRUE) + theme(axis.text.x = element_text(angle = 90, hjust = 1));
print(p);
dev.off();

