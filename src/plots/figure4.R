library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(tidyverse)
library(viridis)
library(forcats)
data <- read.table("data.csv", header=TRUE, sep="\t")
pdf("Figure4-PP_Gene_count_gt1_by_Genome.pdf", width=11.69, height=8.27, pointsize=8, );
data <- data %>%
  gather(key="text", value="value") %>%
  mutate(text = gsub("\\.", " ",text)) %>%
  mutate(value = log(as.numeric(value)));

p <- data %>%
  mutate(text = fct_reorder(text, value)) %>%
  ggplot( aes(x=value, color=text, fill=text)) + geom_histogram(alpha=0.6, binwidth = 0.25) + scale_fill_viridis(discrete=TRUE) +
    scale_color_viridis(discrete=TRUE) + theme( legend.position="none", panel.spacing = unit(0.1, "lines"), strip.text.x = element_text(size = 8) ) + xlab("") +
    ylab("Phylogenetic profile Gene Count > 1 by Genome") +
    facet_wrap(~text);
print(p);
dev.off();

