count_PE_fragments
==================

Count number of paired-end RNAseq fragments mapped to each gene in bacterial genomes
This script is intendent to count the number of sequencing fragments that
overlap each gene (or feature) in the genome. This was written with a bacterial
gene expression pattern in mind meaning there are long transcription units
that contain more than one gene. Since some genes could be small they might be
contained inside operons and thus not be represented in the reads (the ends
of the fragment).
Count the fragment for a gene only if the overlap is higher than a predefined
variable (l). If the fragment is mapped to different locations count each of
these locations.
Input is a gtf (gff) file, read only the lines that are "exon" (column 2),
identifier is quoted after gene_id in column 8.

22/10/14: Added single-end support
07/12/14: Added -r option for reverse reads (useful mainly in single-end)


Usage
=====
  count_PE_fragments.py [options]

Options
=======
--gtf=GTF, -g GTF       GTF file containing the features.
--samfiles=SAMFILES, -s SAMFILES
                        Input sam or bam files, comma separated.
--reverse, -r           Count the genes on the reverse strand of the mapping.
--feature=FEATURE, -f FEATURE
                        Name of features to count on the GTF file (column 2).
--identifier=IDENTIFIER, -i IDENTIFIER
                        Name of identifier to print (in column 8 of the GTF
                        file).
--overlap=OVERLAP, -o OVERLAP
                        Minimal required overlap between the fragment and the
                        feature.
--help, -h              Show this help message and exit.
