# GENPPI parameters
We present the GENPPI possible parameters as a sequential numbered list. The idea is to facilitate the reference of 29 parameters of our software. There is only one mandatory parameter, the directory containing the multifasta files of amino acids. All the optional parameters have their default values. However, one can change as wish.

## Mandatory parameter

### Directory parameter:
1. *-dir* workingdir directory

## Optional parameters

## Parameters of ppi by conserved neighbourhood:

2. *-cnp* conserved neighbourhood score percentage (65-default)
3. *-expt* type of expansion of the conserved gene neighborhood 'fixed' or 'dynamic' (fixed-default)

###  Parameters of conserved neighbourhood by fixed expansion:

4. *-w1* window size 1 (10-default)
5. *-cw1* gene conservation required for window 1 (4-default)
6. *-w2* window size 2 (7-default)
7. *-cw2* gene conservation required for window 2 (3-default)
8. *-w3* window size 3 (5-default)
9. *-cw3* gene conservation required for window 3 (2-default)
10. *-w4* window size 4 (3-default)
11. *-cw4* gene conservation required for window 4 (1-default)

### Parameters of conserved neighbourhood by dynamic expansion:

12. *-ws* dynamic expansion window size (1-default)

## Parameters of ppi by phylogenetic profile:

13. *-ppp* phylogenetic profiles score percentage (30-default)
14. *-ppdifftolerated* difference in phylogenetic profiles tolerated to infer ppi (0-default)

### Amino acid histogram parameter settings for the phylogenetic profile:

15. *-pphistofilter* build the phylogenetic profile of genes with a higher percentage of identity
16. *-ppaadifflimit* amino acid difference limit (0-default)
17. *-ppaacheckminlimit* minimum amount of amino acids to check (26-default)

### Methods of ppi prediction by phylogenetic profile and its parameters:

#### Method 1 - Predict ppi by phylogenetic profile only for interactions predicted by conserved neighborhood

18. *-ppcn*

#### Method 2 - PPi prediction by phylogenetic profile without filters

19. *-ppcomplete*

#### Method 3 - Prediction of ppi by phylogenetic profile with a limit of interactions

20. *-ppiterlimit* maximum number of interactions desired (500000-default)

#### Method 4 - Prediction of ppi by phylogenetic profile with interactions limit by weight

21. *-trim* maximum number of interactions by weight (45000-default)

#### Method 5 - Prediction of ppi by the phylogenetic profile only for genes with profiles that cover a greater or lesser number of genomes than an informed threshold

22. *-threshold* phylogenetic profiles threshold
23. *-plusminus* parameter that receives the greater than or less than sign to apply the ppthreshod filter '' or '' (signs greater than and less than, must be enclosed in single or double quotes)

#### Method 6 - Delete groups of ppi predicted by phylogenetic profile that exceed a limit of interactions by weight

24. *-grouplimit* limit of tolerated interactions to maintain a group of ppi with the same weight (45000-default)

#### Method 7 - To exclude genes with unwanted profiles in predicting ppi by phylogenetic profile

25. *-profiles* number of genomes in the unwanted profiles Entry example: 7 (genes that co-occur in a total of 7 genomes will be excluded). To insert more than one profile, the entry must be enclosed in single or double quotes, and the values separated by semicolons. Example: "7; 15; 21"

## Parameters of ppi by gene fusion:

26. *-genefusion* make ppi predictions by gene fusion
27. *-gfp* gene fusion score percentage (5-default)
Note: if gene fusion is included, there will be an increase of more than 100% in the execution time

## Amino acid histogram parameters:

28. *-aadifflimit* amino acid difference limit (1-default)
29. *-aacheckminlimit* minimum amount of amino acids to check (25-default)

### Possible configurations of the histogram parameters mimicking the Needleman-Wunsch algorithm :

|*-aadifflimit*| *-aacheckminlimit* |minimal protein identity%|
|  :----:      |        :----:      |           :----:        |
|     0        |       26           |             100%        |
|     0        |       25           |             100%        |
|     0        |       24           |            97,96%       |
|     0        |       23           |            96,94%       |
|     0        |       22           |            96,94%       |
|     0        |       21           |            94,68%       |
|     0        |       20           |            91,75%       |
|     0        |       19           |            85,57%       |
|     0        |       18           |            50,00%       |
|     1        |       26           |            97,87%       |
|     1        |       25           |            92,55%       |


