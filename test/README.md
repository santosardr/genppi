##  Command-line Example
To test GENPPI, we provided folders containing two genomes: *Buchnera aphidicola* and *Corynebacterium pseudotuberculosis*. For *Buchenera* genomes, we have an instant test; It should execute speedily. The *Corynebacterium* genomes take a little bit longer to run. However, using the parameter "-expt fixed" for the conserved neighborhood algorithm,  we expect no more than one hour to finish the 50 genome analyzes.

Let's through a step-by-step using the *Buchenera* genomes:

1) Clone the genppi git repository to obtain the entire project on your computer. Having the git installed and using a command-line terminal, you can type:

    `$ git clone https://github.com/santosardr/genppi.git`
2) From a command-line terminal, change the cursor to the folder:

    `$ cd genppi/test/Buchnera_aphidicola`
3) There are two folders inside Buchnera_aphidicola, the *assemblies* and the *refer*. The *assemblies* contains the protein multifasta files for five genomes. A grep looking for the greater than (>) signal tell us the number of proteins per genome:

    `Buchnera_aphidicola$ grep -c '>' assemblies/*.faa`
```text
        assemblies/Ba_Ak.faa:625
        assemblies/Ba_Bp.faa:553
        assemblies/Ba_G002.faa:621
        assemblies/Ba_Sg.faa:620
        assemblies/Ba_Ua.faa:591
```
    However, we recommend working with the *refer* folder. The *refer* contains links for these genomes instead of a raw copy for each file. The purpose is to avoid creating multiple copies of all genomes because GENNPI writes his results in the form of subfolders and files within the folder passed as a parameter at runtime. If you execute GENPPI twice in the same folder, GENPPI will overwrite previous results. Unfortunately, for MS Windows users, the *refer* folder is useless because there is no way to create links like Mac and Linux do; copy the assemblies folder or use the current one (GENPPI do not alter the fasta files).

4) Make a copy of the data source folder. Here the destination is the test1 folder.
For Linux and Mac: 

    `$ cp -r refer test1`

    The "-r" should guarantee link preservation during copies, saving your disk space. Pay attention to the fact that offered links are relative to one folder up.
    For MS-Windows:

    `$ cp -r assemblies test1`

5) If the folder copying was successful, now it's time to execute GENPPI. 

    5.1) *Buchenera* genomes are one of the smallest genomes we have ever known; we expect a fast GENPPI execution and a few hundred edges for each genome. 

    `$ genppi -dir test1/`

    Checking the numerical results:

    `$ wc -l test?/ppi-files/*.sif`
```text
          435 test1/ppi-files/Ba_Ak.sif
           90 test1/ppi-files/Ba_Bp.sif
         **462 test1/ppi-files/Ba_G002.sif**
          243 test1/ppi-files/Ba_Sg.sif
          350 test1/ppi-files/Ba_Ua.sif
```
    For standard GENPPI run (default parameters), the Ba_G002 has the more extensive interaction network with 462 edges comprising 126 proteins. Just looking for this output on the screen, we cannot know about the number of unique proteins. Please, hold on to that. Soon I will show you how.

    5.2) I will relax the GENPPI parameters to obtain more interactions in the final networks. When looking for conserved phylogenetic profiles, I meant to do such a relaxing parameter telling the GENNPI algorithms to accept as similar proteins those with at least 50% of identity (-aadifflimit 0 -aacheckminlimit 18). Additionally, I will ask for a dynamic expansion in the conserved neighborhood algorithm; it will start with a minimum window size of four to infer conservation (-ws 4). If the algorithm is successful for an initial ws, it will expand the window size by four units for subsequent well-success expansions. Besides, I will ask GENPPI for not using any filter for phylogenetic profiles (-ppcomplete).  After all, this is our command:

    `$ genppi -expt dynamic -ws 4 -ppcomplete -aadifflimit 0 -aacheckminlimit 18 -dir test2/`

    Checking the numerical results:

    `$ wc -l test?/ppi-files/*.sif`
```text
    435  test1/ppi-files/Ba_Ak.sif
    90   test1/ppi-files/Ba_Bp.sif
    462  test1/ppi-files/Ba_G002.sif
    243  test1/ppi-files/Ba_Sg.sif
    350  test1/ppi-files/Ba_Ua.sif
    1232 test2/ppi-files/Ba_Ak.sif
    413  test2/ppi-files/Ba_Bp.sif
    1263 test2/ppi-files/Ba_G002.sif
    821  test2/ppi-files/Ba_Sg.sif
    1002 test2/ppi-files/Ba_Ua.sif
```

    Our new interaction network for Ba_G002, compared to the one obtained with GENPPI default parameters, has 1263 edges comprising 172 proteins or a three-fold in the number of interactions and a 36% increase in the number of unique proteins.

## Exploring GENPPI results

Once we have some networks created by GENPPI, what can we further do with those? Well, this answer will depend on your researching needs. However, the most common requirement is to visualize the results.  Our purpose was to produce textual outputs easy to deal with and primarily used. That's the case for the DOT format, broadly used in the GEPHI tool, SIF extensively used by Cytoscape, and R plugins capable of reading Cytoscape and other tabular formats. **GENPPI outputs two data formats,  DOT and SIF, for each genome (multifasta file) analyzed**. In the case of R as the visualizing/processing tool for interaction networks, one should remember that the SIF format is a three-column file where the first and third columns are the interacting protein identifications. So, if one needs a Tab Separated Value (TSV) file or Comma Separated Value (CSV) file, import the SIF file in a spreadsheet program and export only the first and third columns of the SIF file. Another possibility is to run a command line like this one:

    $ cut -f 1,3 test2/ppi-files/Ba_G002.sif
    "BUMPG002_CDS00574"	"BUMPG002_CDS00573"
    "BUMPG002_CDS00574"	"BUMPG002_CDS00572" (and so on)

Or, to redirect the TSV for a file, just type:

    $ cut -f 1,3 test2/ppi-files/Ba_G002.sif > Ba_G002.tsv

We will give you basic instructions to load the interaction networks created by GENPPI in some visual tools. Please, check the appropriate software documentation about how to install each software.

### GEPHI https://gephi.org/

A fast way to open a DOT file with GEPHI is by calling it in the command line. 

    $ gephi test2/ppi-files/Ba_G002.dot &

The GEPHI interface will ask you some basic questions about open a new workspace or append the data to an existing one. It also will ask you about the edges merge strategy. I use to sum several edges of the same pair of interacting proteins. For GENPPI, which could create three edges for a couple of interacting proteins, GEPHI will summarize in only one connection arrow. The OK button will open the interaction network. Initially, the interaction network appearance is a mass, completely aleatory distribution of vertices and connections. A more elegant view, for instance, is obtained by the options Window->Layout->Yifan Hu. One can see clusters of vertices. Another possibility is to calculate statistics according to the current topology. The Window->Statistics open a lateral menu of possible calculations. For instance, the Network Diameter option allows us to calculate the Betweenness Centrality and other measures. All the calculated data is available in a spreadsheet-like format in the Data Laboratory view.

![GEPHI sample](https://github.com/santosardr/genppi/raw/master/doc/Ba_G002-GEPHI.png)

### Cytoscape https://cytoscape.org/

After opening Cytoscape, generally via mouse action, you should ask for File->Import->Networking from File.
Navigate to the folder used for GENPPI output and load the file "Ba_G002.sif". The Style menu (located at the left lateral menu) allows you to customize your view. To make a graph similar to the GEPHI style, you can change the Shape and Width options to Ellipse and 30, respectively. The Tool->Analyse Network menu allows for several topology measures at once, including Betweenness Centrality. Depending on the number of edges and vertices in your network, running this option can be very time-consuming.

![Cytoscape sample](https://github.com/santosardr/genppi/raw/master/doc/Ba_G002-Cytoscape.png)

### R https://www.r-project.org/

Install R packages to handle protein interaction networks sometimes require previous dependencies. For instance, On installing the package *graph* demands the install of the installer package of the software provider, the *BiocManager*. This process takes several minutes of downloading and installing:


    $ sudo R
```R
    if (!requireNamespace("BiocManager", quietly = TRUE)) \
    install.packages("BiocManager")
    BiocManager::install("graph", version = "3.8")
    q()
```

 After finishing the installation process, I will use chapter 11 of "A Little Book of R for Bioinformatics!" written by  Avril Coghlan from the Trust Sanger Institute, Cambridge, to demonstrate how to use GENPPI output in R. 
The main trick to compute interaction networks by GENPPI in R statistical software is creating a data frame with GENPPI output. Fortunately, Mr. Coghlan provided us an [**R function**](https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter11.html#reading-in-protein-protein-interaction-data-in-r) (click the bold text to navigate) to read a two-column file as an R data frame. The only modification I had to perform here to load my data was at line 8. Instead:
```R
    protnames <- c(levels(proteins1),levels(proteins2))
```
I need to change it to:
```R
    protnames <- c(proteins1,proteins2)
```

Now, we can follow most of the examples of "A Little Book of R for Bioinformatics!". Let us start by reading the data. In Linux and Mac, we will create a TSV from a SIF file generated by GENPPI:

    $ cut -f 1,3 test1/ppi-files/Ba_Sg.sif > Ba_Sg.tsv

Please, pay attention that now I'm using a smaller network created when I ran the default parameters of GENPPI. Now let us open the R statistical software in the same folder where you saved the TSV:

    $ R

Please copy and paste the *makeproteingraph* function in a text editor and update line 8, as I pointed above. After the edit, copy and paste the edited function in the R console. Now, load the data file:
```R
    BaSg <- makeproteingraph("Ba_Sg.tsv")
```
After you navigate to the top of the web page of "A Little Book of R for Bioinformatics!" you will find several cool things to make with our network. For instance, this function shows the adjacent nodes to one specific:
```R
    adj(BaSg,"BUsg_616")
    $BUsg_616
    [1] "BUsg_613" "BUsg_573" "BUsg_574" "BUsg_575" "BUsg_576"
```
We can show the degree distribution of the graph:

```R
    mydegrees <- graph::degree(BaSg)
    sort(mydegrees)
    BUsg_018 BUsg_590 BUsg_595 BUsg_019 BUsg_065 BUsg_020 BUsg_591   BUsg20 
    1        1        1        2        2        2        2        2 
    (and so on)
    BUsg_393 BUsg_394 BUsg_395 BUsg_396 BUsg_397 BUsg_607 BUsg_470 BUsg_471
    10       10       10       10       10       10       10       10 
    BUsg_473 BUsg_613 BUsg_601 BUsg_604 BUsg_603 
    10       10       11       12       13 
```
We can plot a histogram of the degree distribution:

```R
    hist(mydegrees, col="red")
```
![R sample 1](https://github.com/santosardr/genppi/raw/master/doc/BaSg-R1.png)

To Finish, we can plot the network interaction on the screen:

```R
    library("Rgraphviz")
    mygraphplot <- layoutGraph(BaSg, layoutType="fdp")
    renderGraph(mygraphplot)
```

![R sample 2](https://github.com/santosardr/genppi/raw/master/doc/BaSg-R2.png)

Thank you for your patients. 

Enjoy it.
