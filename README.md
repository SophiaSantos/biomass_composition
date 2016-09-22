## Biomass Composition

The java application created is a tool that provides the estimation of biomass composition in nucleotides and amino acids, with input files containing sequences from DNA, RNA and protein, in the FASTA format. When expression data are available, it can also be used, provided in a csv file containing percentages of each gene/protein. The output of this tool are the amino acid, nucleotide and deoxynucleotide compositions in percentage and in mmol/gDW. These can be directly included in the biomass equation.

To obtain the results it is only necessary to click in the “Determine” button. It is also possible to export the obtained data to a file in csv format, by clicking in the “Export” button.
 
All data obtained can be easily exported to a csv file.

This application allows to obtain the results rapidly and is also a user-friendly tool for users with any or little background in informatics.

### Obligatory inputs

- Input files with sequences of Proteins, DNA and RNA, exclusive in the FASTA format; 

- Transcriptomic data, if available, in csv format, with two columns separated by semicolon: the first column should contain gene identifiers and the second the expression factor in percentage. In this case, the FASTA file with protein sequences should have the same gene identifiers at the beginning of the sequence header;

- Percentage (number between 0 and 1) of each type of RNA (mRNA, rRNA and tRNA);

- Percentage of the cellular content in each macromolecule (Protein, DNA and RNA) in percentage (number between 0 and 1).


### Prerequisites

[Git][] and [JDK 8 update 20 or later][JDK8 build]

[Git]: http://help.github.com/set-up-git-redirect
[JDK8 build]: http://www.oracle.com/technetwork/java/javase/downloads