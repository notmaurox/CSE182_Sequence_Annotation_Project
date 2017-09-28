# CSE182_Sequence_Annotation_Project
Some deliverables from a sequence annotation project I did for my Biological Databases class where we analyzed a portion of Acinetobacter baumanni’s proteome

# annotation.py
This python script takes in a .fasta file of query sequences, in the context of this project we used protein sequences found in Acinetobacter baumanni, and searches them against 3 different data bases. It filters the results and pulls out interesting information from the best matches to each database and writes output files. In this project we formatted our output to fit the needs of other project groups in the class building online web interfaces to view our data. Our data was shared with 4 other groups so our script had to write 4 files each matching the specifications of each group’s web platform. 

Databases and tools used:
BLASTp on UniProt(Reviewed Swiss-Prot)Data Base / UniProt database downloaded and BLASTp run locally
HMMScan on Pfam
ExPASy ScanProsite package

# fileUpdater.py
This script was used to update all the data files that would go to the other project teams by incorporating our hand written comments on each annotation. 

# harambaesAnnotationsFinal.tsv
An example of the output we gathered that was provided to other project teams.

# rawProtienAnnotationDataIndexd.json
As an extra credit option, we wrote all our raw search results and match data to an indexed json file.

# CSE182FinalProjectReport.pdf
A complete description of the workflow used in the project, implementation details, interesting findings, and a reflection on the process as a whole.

# sequences.fasta
File containing query sequences provided by Dr Vineet Bafna in CSE 182 during Spring of 2017
