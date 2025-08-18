import os
import csv
from glob import glob
from Bio import SeqIO
from fuzzysearch import find_near_matches

def find_cys2_motif(sequence, similarity_threshold, cys2_motif="ccaagcacctgccgtca"):
    sequence_lower = sequence.lower()
    cys2_motif_lower = cys2_motif.lower()

    matches = find_near_matches(cys2_motif_lower, sequence_lower, max_l_dist=int(len(cys2_motif) * (1 - similarity_threshold)))
    filtered_matches = [match.start for match in matches]
    return filtered_matches

def find_smacla_motif(sequence, similarity_threshold):
    sequence_lower = sequence.lower()
    smacla_motif = "cccgggcggcttgtcttttaaaggtttgcaaggcgggcggggtcgtccgttccggtggaaataatatatcgat"
    smacla_motif_lower = smacla_motif.lower()

    match = find_near_matches(smacla_motif_lower, sequence_lower, max_l_dist=int(len(smacla_motif) * (1 - similarity_threshold)))
    
    if match and match[0].dist <= int(len(smacla_motif) * (1 - similarity_threshold)):
        return match[0].end
    return -1

def find_stop_codon(sequence, start_index, search_range=250):
    stop_codons = ["taa", "tag", "tga"]
    search_start = start_index
    search_end = min(start_index + search_range, len(sequence))
    search_window = sequence[search_start:search_end].lower()
    for stop_codon in stop_codons:
        stop_codon_pos = search_window.find(stop_codon)
        if stop_codon_pos != -1:
            return start_index + stop_codon_pos
    return -1

def generate_copy_files_from_sequence(record, similarity_threshold, output_folder, csv_file):
    sequence = str(record.seq)
    cys2_motif = "ccaagcacctgccgtca"
    
    cys2_start_indices = find_cys2_motif(sequence, similarity_threshold, cys2_motif)

    print(f"For record {record.id}: Number of Cys2 motifs found: {len(cys2_start_indices)}")
    print(f"For record {record.id}: Cys2 start indices: {cys2_start_indices}")

    num_fragments = 0
    num_smacla = 0
    prev_stop_codon = -1

    for copy_num, start_index_cys2 in enumerate(cys2_start_indices, start=1):
        start_index_smacla = find_smacla_motif(sequence, similarity_threshold)
        start_index_stop_codon = find_stop_codon(sequence, start_index_cys2 + len(cys2_motif))

        if start_index_smacla != -1:
            num_smacla += 1
            output_file_name = os.path.join(output_folder, f"{record.id}_copy{copy_num}.fasta")
            with open(output_file_name, 'w') as output_file:
                output_file.write(f">{record.id}_copy{copy_num}\n")
                output_file.write(sequence[prev_stop_codon + 1:start_index_smacla])  
            num_fragments += 1
            prev_stop_codon = start_index_smacla
        elif start_index_stop_codon != -1:
            num_fragments += 1
            output_file_name = os.path.join(output_folder, f"{record.id}_copy{copy_num}.fasta")
            with open(output_file_name, 'w') as output_file:
                output_file.write(f">{record.id}_copy{copy_num}\n")
                output_file.write(sequence[prev_stop_codon + 1:start_index_stop_codon + 3])  
            prev_stop_codon = start_index_stop_codon

        update_csv(csv_file, record.id, start_index_cys2, num_fragments)

    print(f"For record {record.id}: Number of identified smacla motifs: {num_smacla}")
    print(f"For record {record.id}: Number of generated fragments: {num_fragments}")

def update_csv(csv_file, sequence_id, relative_position, num_copies):
    existing_data = read_existing_csv(csv_file)


    if sequence_id in existing_data:
        existing_data[sequence_id]["RelativePosition"] = relative_position
        existing_data[sequence_id]["NumCopies"] += num_copies
    else:
        existing_data[sequence_id] = {"RelativePosition": relative_position, "NumCopies": num_copies}


    with open(csv_file, 'w', newline='') as csvfile:
        fieldnames = ['SequenceID', 'RelativePosition', 'NumCopies']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for sequence_id, data in existing_data.items():
            writer.writerow({'SequenceID': sequence_id, 'RelativePosition': data['RelativePosition'], 'NumCopies': data['NumCopies']})

def read_existing_csv(csv_file):
    existing_data = {}
    try:
        with open(csv_file, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                sequence_id = row['SequenceID']
                relative_position = int(row['RelativePosition'])
                num_copies = int(row['NumCopies'])
                existing_data[sequence_id] = {'RelativePosition': relative_position, 'NumCopies': num_copies}
    except FileNotFoundError:
        pass

    return existing_data

if __name__ == "__main__":
    input_pattern = "pilS"  
    output_folder = "."  
    similarity_threshold = 0.76
    csv_file = "pilS_report.csv"  


    input_files = glob(f"*{input_pattern.lower()}*")

    for input_file in input_files:
        print(f"Processing input file: {input_file}")
        records = list(SeqIO.parse(input_file, 'genbank'))

        for record in records:
            print(f"Processing record: {record.id}")


            generate_copy_files_from_sequence(record, similarity_threshold, output_folder, csv_file)
