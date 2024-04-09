install.packages("BiocManager");
BiocManager::install("graph");
library("graph")

makeproteingraph <- function(myfile)
  {
     # Function to make a graph based on protein-protein interaction data in an input file
     require("graph")
     mytable <- read.table(file(myfile)) # Store the data in a data frame
     proteins1 <- mytable$V1
     proteins2 <- mytable$V2
	protnames <- c(proteins1,proteins2)
     # Find out how many pairs of proteins there are
     numpairs <- length(proteins1)
     # Find the unique protein names:
     uniquenames <-  unique(protnames)
     # Make a graph for these proteins with no edges:
     mygraph <- new("graphNEL", nodes = uniquenames)
     # Add edges to the graph:
     # See http://rss.acs.unt.edu/Rdoc/library/graph/doc/graph.pdf for more examples
     weights <- rep(1,numpairs)
     mygraph2 <- addEdge(as.vector(proteins1),as.vector(proteins2),mygraph,weights)
     return(mygraph2)
  }
system("cut -f 1,3 Bvmandacaium.sif > Bvmandacaium.tsv")
bvm <- makeproteingraph("Bvmandacaium.tsv")
mydegrees <- graph::degree(bvm)
summary(mydegrees)
bvm
