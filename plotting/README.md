# Graphical GENPPI's outputs

We opt not to try to develop code for drawing bar and box plots. The reason is 
simple: tons of excellent libraries in R, Python, Matlab, and others. It allows 
us to focus on our best: producing reliable interaction networks.
This is an example of a possible graphical plotting you can do from GENPPI's output:

![ Phylogenetic profile Box plot for Corynebacterium genera](https://github.com/santosardr/genppi/raw/master/plotting/sample.pdf)

**NOTE**
If you are not interested in continuing this tutorial but do want to see plots 
from your data, we also offer you a web site for this task.

**<http://bioinformatics.college/genppivisual/>**

All you need is to upload the files named "report.txt" from the data created by 
GENPPI. The limitation is the size of the files. 

The main script you should execute to create the plots is the *genppivisual*. Its 
execution is quite simple since their R's dependencies are correctly installed.
We have five R's scripts called by our main bash script *genppivisual*, one for 
each plot derivate from the GENPPI reports. To a success report plot, we must a 
few R libraries previously installed. Despite being a few, some of these 
libraries depend on several others. Because of that, one should expect several 
(~30) minutes to finish all the dependencies set, plus at least 2 GB of RAM.
Here are the steps for installation and execution:

1 - R statistical software and development packages.
For the Linux Ubuntu OS, one command can do this job.

> $ sudo apt install *package*

However, you need to install each software separately due to dependencies requirements, where *package* is:

+ r-base
+ libxml2-dev
+ libfontconfig1-dev
+ libcurl4-openssl-dev
+ libssl-dev
+ libcairo2-dev
+ libgit2-dev

2 - Libraries.

$ sudo R

> \>install.packages(c("dplyr", "forcats", "ggplot2", "hrbrthemes", "tidyr",
"tidyverse", "viridis"), dependencies=TRUE)

Please, check if each library has been corrected installed. 
Try to load them separately.

3 - Exit from R

> \>q()

4 - awk software

I use awk to process the reports and format appropriate files for R.

> $ sudo apt install gawk

4 - Finally, we can plot a bar and box plots. Each folder created by GENPPI has 
a "report.txt" file, which we can use to plot. Just indicate one of the files 
named "report.txt". For instance, switch to the "visual" folder and type:

> $ ./genppivisual $PWD/../../test/Corynebacterium_pseudotuberculosis/sample-report/phylogenetic-profiles-report/report.txt

A word about this command: the $PWD ensures to find the relative path two 
folders up. Otherwise, we could receive a *"file not found"* error. For your use, 
I recommend passing the full file path.
After successful execution, you will find two PDF files in a folder named 
"report.txt.tmp" in the same folder.

We listed the R's scripts and respective producing actions.    

* GENPPI's folder --> phylogenetic-profiles-report

| script    | output                                  |
|-----------|-----------------------------------------|
| figure1.R | Figure1-BoxPlot-PP-by-Genomes.pdf       |
| figure4.R | Figure4-PP_Gene_count_gt1_by_Genome.pdf |

**Pay attention**
GENPPI only creates the phylogenetic-profiles-report folder if you use the  
-ppcomplete parameter.

* GENPPI's folder --> gene-neighborhood-conservation-report

| script    | output                                  |
|-----------|-----------------------------------------|
| figure3.R | Figure3-BoxPlot-CN-by-Genomes.pdf       |

* GENPPI's folder --> ppi-report

| script    | output         |
|-----------|----------------|
| figure2.R | Figure2-CN.pdf |
| figure5.R | Figure5-PP.pdf |

Enjoy it.

Anderson Santos
