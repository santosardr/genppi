# genppi
Ab initio protein interation network generator

GENPPI is a software written in Common Lisp and compiled with SBCL to create 
complex networks of predicted proteins from tens or hundreds of genomes. Your 
main limitation is the amount of RAM in your server. For instance, using a 
conventional computer of 8 Gigabytes of RAM, the software can deal with at least
50 genomes, each possessing on average 2200 proteins (Corynebacterium). However,
this same machine cannot process 80 genomes with an average of five thousand 
proteins (Escherichia coli). To achieve this last task well done, we need to 
compile GENPPI to use at least 8 Gigabytes of RAM and use the program version 
capable of storing the bulk data on disk instead of memory, the GENPPIDB executable.
For now, I had compiled all versions available in the binaries folder to use 
different sizes of RAM. While I do not release the source code, if you need a 
version compiled with more RAM, please get in touch with me. SBCL includes the 
whole core of libraries for each executable, and it explains the size of the 
software in megabytes. 
Till today, I haven't a graphical interface for GENPPI since I conceived it as 
a command-line tool for the Linux OS. Calling the software without any arguments
results in it printing all the possible parameters and their combinations. It 
also happens when starting the program with the -help option. 
I will show a fast track to obtain results with this software. First of all, 
(i) create a folder containing multi-fasta files of predicted proteins for your 
genomes. We will treat each file as a unique genome. I recommend you name these 
files significantly because a lot of reports will mention these names. I also 
suggest keeping the protein names as simple as possible and short. I have 
another program that you can use for this task: I called it valifasta. Valifasta
tries to figure out a combination of strings identifying without redundancy all
protein within a multi-fasta file. If such an assignment does not succeed, 
valifasta creates a numerical and sequential identification for each protein and
uses it instead of the original proteins.
To run the program (ii), there are a lot of parameters and combinations of those. 
I recommend you to start simple. Probably you will not reach out to a 
comprehensive interaction network on your first try. As you master the simple, 
do further steps. As a rule of thumb, a good interaction network should 
represent the majority of your proteins (number of vertices), possess at least 
two or three thousand interactions (number of edges), and on average a few 
hundred (<200) interactions per protein. Some proteins can have thousands of 
connections, but these should not be in a representative number. A configuration
that I use to frequently is this one:
genppidb -ppcomplete -expt fixed -w1 7 -cw1 4 -ppdifftolerated 1  -pphistofilter
-dir  summary/

Meaning:
"-ppcomplete" is an excellent parameter to run at least once for each set of 
genomes. It does not restrict the number of interactions concerning  
Phylogenetic Profile (PP); it just let it go. Depending on a small number of 
genomes in your folder, the consequence is a massive and undesirable number of 
edges as a final result. The reason is that for a few genomes, very closely 
related, the majority of the phylogenetic profiles will be very conserved. 
Try to analyze related genomes but not necessarily very similar to diminish the 
number of possible expected interactions. The phylogenetic profile report only 
is created when we settle this parameter. If you decide to restrict the number 
of connections, do not use this parameter but limit the interactions. Optional 
parameters to limit edges are -ppiterlimit, -trim, and -threshold, just citing 
some of them. Check the -help options for a full list.

"-expt fixed" set the program for comparing neighbors genes at most "-w1" genes,
and using "-cw1" as the smaller number of genes for concluding a Conserved 
Neighbourhood (CN). CN and PP are the primary methods to map significant 
connections in an interaction network. The counterpart of "-expt fixed" is 
"-expt dynamic -ws 3", for instance. It conducts a systematic expansion of a 
windows' limits on analyzing a CN. However, dynamic growth has time-consuming 
proportional to -ws value.

"-ppdifftolerated" give you more interactions as a final result of the 
Phylogenetic Profile (PP) analyses. This parameter's integer value defines the 
level of tolerance to accept two PPs as similar between genomes.

"-pphistofilter" set the GENPPI to be too restrictive when deciding the 
similarity between two proteins, but only during the PP analyses. The lack of 
this parameter can also create a much larger number of edges at the final 
interaction network.

"-dir" is the location and name of the folder where you deposited the 
multi-fasta files to be processed.

Enjoy it.

Anderson Santos
santosardr@ufu.br