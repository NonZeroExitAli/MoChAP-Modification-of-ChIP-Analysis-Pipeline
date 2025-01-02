from Local_alignment import *
def generate_sam_sw(header1, aligned_seq1, header2, reference_start, is_reverse=False, mapping_quality=40):
    bitwise_flag = 16 if is_reverse else 0  # Set the flag to 16 if the read is reverse complemented

    # Write the sam entry in a shape of SAM file
    sam_entry = ""
    sam_entry += f"{header1.split()[0]}\t{bitwise_flag}\t{header2.split()[0]}\t1\t{reference_start}\t{mapping_quality}\t{len(aligned_seq1)}M\t*\t0\t0\n{aligned_seq1}\t*\n"

    return sam_entry
# FASTA & FASTQ files
Sample_file = "MG1655-RpoN-1_R1s.fastq"
Reference_file = "NC_000913-complete_genome.fasta"

headers1, sequences1 = read_fasta(Sample_file)
headers2, sequences2 = read_fasta(Reference_file)

# Access the sequence from the reference file
reference_sequence = sequences2[0]

# Initialize the SAM file main header
sam_entry = "@HD VN:1.6 SO:coordinate\n"
# Align all sequences from sample file( to the reference sequence
for i, (header1, sequence1) in enumerate(zip(headers1, sequences1)):
    aligned_seq1, aligned_seq2, reference_start = smith_waterman(sequence1, reference_sequence)

    # Add the data as SAM entry
    sam_entry += generate_sam_sw(header1, aligned_seq1, headers2[0], reference_start)

    with open('output.sam', "w") as sam_file:
        sam_file.write(sam_entry)