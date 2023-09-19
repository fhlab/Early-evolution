# Early-evolution

Custom python scripts used for data analysis in https://www.biorxiv.org/content/10.1101/2023.02.13.528392.abstract

Both scripts take a .tsv file from the DiMSum pipeline (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02091-3) as input and analyse sequences and read counts. 

## Deletion Histogram Generator

The deletion histogram generator creates a histogram and .csv file of the occurance of single nucleotide deletions with a variant merge .tsv file from the DiMSum pipeline and either "r1"  or "r2" as round ID as input:

'python histogram_generator_deletions.py --variant_merge_file VARIANT_MERGE_FILE --round ROUND_ID'

## Stop Histrogram Generator

The stop histogram generator generates a histogram for the occurence of stop codons and for a residue of interest in the protein sequence. It takes a variant merge .tsv file from the DiMSum pipeline, a round ID and a residue of interest as input: 

'python stop_and_residue_of_interest_histogram_generator.py --variant_merge_file VARIANT_MERGE_FILE --round ROUND_ID --reside RESIDUE_1_LETTER_CODE'
