library(tidyr)
library(ggplot2)
data=read.csv("data.csv",header=TRUE,stringsAsFactors=FALSE, colClasses =c("character", "integer", "integer"));
pdf("Figure5-PP.pdf", width=11.69, height=8.27, pointsize=8, );
p<-ggplot(data, aes(x=Genome, y=PP ))+ geom_bar(stat = "identity", position = 'dodge', fill="#66FAE9") + theme(axis.text.x = element_text(angle = 90, hjust = 1));
print(p);
dev.off();

