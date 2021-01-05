library(ggplot2)
data=read.csv("cn.csv",header=TRUE,stringsAsFactors=FALSE, colClasses =c("character", "integer"), sep="\t");
pdf("Figure3-BoxPlot-CN-by-Genomes.pdf", width=11.69, height=8.27, pointsize=8, );
p<-ggplot(data,aes(x=Genome, y=Genes))+geom_boxplot(varwidth = TRUE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ labs(x="Genome", y="Gene count by CN");
print(p);
dev.off();

