library(tidyr)
library(ggplot2)
data=read.csv("data.csv",header=TRUE,stringsAsFactors=FALSE, colClasses =c("character", "integer", "integer"));
pdf("Figure2-CN.pdf", width=11.69, height=8.27, pointsize=8, );
p<-ggplot(data, aes(x=Genome, y=CN ))+ geom_bar(stat = "identity", position = 'dodge', fill="#FF6666") + theme(axis.text.x = element_text(angle = 90, hjust = 1));
print(p);
dev.off();

