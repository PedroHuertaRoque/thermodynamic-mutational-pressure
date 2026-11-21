# Thermodynamic mutational pressure

Analysis scripts for a study of GC-changing substitutions during serial *in vitro* DNA amplification. The study analyzed 62 Nanopore-sequenced samples from nine amplification series involving three targets (*evgA*, *ygjV*, and *phnH*) and three amplification systems (DreamTaq, the Klenow fragment, and recombinase polymerase amplification, RPA).

## Analysis workflow

| Step | Script | Analysis stage |
| --- | --- | --- |
| 1 | `1.Mapping.sh` | Read alignment against the reference sequence. |
| 2 | `2.Filter.sh` | Filtering of alignments for downstream analysis. |
| 3–4 | `3.Coverage1.sh`, `4.Coverage2.sh` | Sequencing coverage and depth summaries. |
| 5 | `5.Count.sh` | Nucleotide counting. |
| 6 | `6.Selection.sh` | Selection of candidate substitutions. |
| 7 | `7.Mutations.py` | Evaluation of candidate persistence across amplification rounds. |
| 8 | `8.Allele.py` | Allele-frequency analysis of retained substitutions. |
| 9 | `9.GC_percentage.py` | Calculation of the GC index from retained substitutions. |

## Inputs and software

The workflow starts from demultiplexed Nanopore FASTQ files. It also uses the study reference sequences and the mapping between sample identifiers, amplification systems, target genes, and rounds. Check the input paths and settings at the top of each script before running it. The workflow uses Bash, Python, Minimap2, and SAMtools. The manuscript's Materials and Methods section describes the filtering and variant-retention criteria.

## Supplementary data

- **Online Resource 1:** Reference amplicon sequences, reaction compositions, sequencing quality and alignment summaries, the RPA candidate catalog (Table S11), background allele-frequency estimates, and supplementary figures.
- **Online Resource 2 (Table S9):** Retained DreamTaq candidate catalog.
- **Online Resource 3 (Table S10):** Retained Klenow candidate catalog.
