# Motif Finder: Identification et extraction of genetic motifs

This script uses two repeats are characteristic of the pilE and pilS loci to allows the identification and the extraction pilS loci. 
It works by identifying  specific genetic motifs (cys2 and  SmaCla) in genomic sequences, and extracts the surrounding sequences fragments pilS and pilE. Then the extracted fragments are stored in a fasta file, and a CSV file report is generated. 
---

input_pattern: word (pattern) that has to be in the name to be analyzed, (by default, "pilS") 
output_folder: Output folder for the extracted FASTA files (by default ., the source rfolder) 
similarity_threshold: Threshold of similarity for the motifs search (by dfaut 0.76) 
csv_file: name of the CSV report by default "pilS_report.csv"

----
Input files: 
Genbank files named with the format pilS_sample1.gb, pilS_sample2.gb, etc.

Output files: 
- FASTA files: For each fragment, a FASTA File is generated with the  format {id_sequence}_copy{numéro}.fasta

- CSV file: a report that contains for each sequence
--SequenceID: Id of the original sequence
--RelativePosition: Relative position of the cys2 motif
--NumCopies: Number of generated fragments 

---

 Dependencies 
- Python 3.x
- Python Libraries:
--`Biopython` 
--`fuzzysearch` 
--`csv` et `glob` (included in standard python library)

To install the missing dependencies: 
pip install biopython fuzzysearch

