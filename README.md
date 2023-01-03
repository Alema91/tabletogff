# Tabletogff
From table annotations create gff.

## Overview
These are some scripts to convert various features and annotations into a GFF-like file for use in genome browsers. There are also general purpose tools for genome annotation.

GFF (generic feature format) is a tab-delimited pseudoformat for annotating genomes, consisting of 8 well-defined columns, and a free-text 9th column (with some "guidelines" of what to include). Most work in fixing a GFF involves fixing the final column. While this is not a great format, most formats are insufficient for one reason or another. Ultimately, it would be better if certain fields were explicit, like ID (a unique machine-understandable identifier, which could be accession, etc.), Parent (any other ID), and Name (like gene name). This would also allow for fast indexing or sorting by ID, Parent, or Name, which would also make it easier to extract or process genes or features for graphics.

## Annotations by columns:

| Column | Name       | Description                                                                                                                                                                                        |
|--------|------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1      | seqid      | The sequence ID, usually this is the contig, scaffold or chromosome (e.g. ctg0123, scaffold_99, chrX), and NOT the name of the gene or feature (which goes in column 9)                            |
| 2      | source     | Usually the program that generated the data (e.g. blastp, AUGUSTUS)                                                                                                                                |
| 3      | type       | Type of feature (e.g. gene, exon), typically using standarized terms                                                                                                                               |
| 4      | start      | First base in the feature                                                                                                                                                                          |
|      5 | end        | Last base the feature                                                                                                                                                                              |
|      6 | score      | Float (of an arbitrary scale), could be something like coverage (say 0-1000), percent identity (0.0-100.0)                                                                                         |
|      7 | strand     | Which DNA strand, as forward (+), reverse (-), unstranded (.) or unknown (?)                                                                                                                       |
|      8 | phase      | Phase of the coding sequence as 0, 1, or 2 (i.e. whether the exon ends mid-codon), only applies to CDS features                                                                                    |
|      9 | attributes | All other information, as a string of pairs of key=value;, though this is variable depending on version or program. Most problems with GFF are problems with parsing information from this column. |

In our test data:

| N | Test data     | Template data | Explanation                                                                                                                                               |
|---|---------------|---------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1 | Gene ID       | seqid         | The sequence ID, usually this is the contig, scaffold or chromosome (e.g. Chromosome I), and NOT the name of the gene or feature (which goes in column 9) |
| 2 | Scaffold name | type          | Type of feature (e.g. gene, exon), typically using standarized terms                                                                                      |
| 3 | Gene start    | start         | First base in the feature                                                                                                                                 |
| 4 | Gene end      | end           | Last base the feature                                                                                                                                     |
| 5 | ORF Lenght    | length        | length feature                                                                                                                                            |
| 6 | ORF start     | start         | First base in aac                                                                                                                                         |
| 7 | ORF end       | end           | Last base in aac                                                                                                                                          |
| 8 | Strand        | strand        | Which DNA strand, as forward (+), reverse (-), unstranded (.) or unknown (?)                                                                              |
| 9 | mRNA sequence | attributes    | sequence; InterProScan; Pfam; InterPro     |