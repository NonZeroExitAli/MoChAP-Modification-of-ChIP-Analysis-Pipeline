import os
import requests

#This is the file where I validate and work on the data (MA0047.3.sites)
# def read_sequences_from_file(file_path):
#     original_sequences = []
#
#     with open(file_path, "r") as file:
#         for line in file:
#             if not line.startswith(">"):
#                 original_sequence = line.strip().upper()
#                 original_sequences.append(original_sequence)
#
#                 if len(original_sequences) == num_sequences:
#                     break
#
#     return original_sequences

#This is the conversion code for the file format that I will obtain from Ali (real data)

def SAM_Read(file_name):
    sequences = []
    with open(file_name, "r") as samf:
        for line in samf:
            if line.startswith("@"):
                continue
            sequences.append(line.strip())
    return sequences

def calculate_frequencies(sequences):
    nucl = set(''.join(sequences))
    total_bases_forward = len(sequences) * len(sequences[0])
    freq_dict_forward = {base: sum(seq.count(base) for seq in sequences) / total_bases_forward for base in nucl}
    return {'forward': freq_dict_forward}

def calculate_motif_matrix(sequences):
    nucl = set(''.join(sequences))
    motif_length = len(sequences[0])
    motif_matrix_forward = np.zeros((motif_length, len(nucl)))

    for i in range(motif_length):
        for j, base in enumerate(nucl):
            motif_matrix_forward[i, j] = sum(seq[i] == base for seq in sequences) / len(sequences)

    return {'forward': motif_matrix_forward}

def conserved_motif_file_exists(motif):
    conserved_file_path = f"conserved_motifs/{motif}_conserved.txt"
    return os.path.exists(conserved_file_path)

def check_conserved_motif(sequences, motif):
    motif_locations = []

    if conserved_motif_file_exists(motif):
        with open(f"conserved_motifs/{motif}_conserved.txt", "r") as motif_file:
            conserved_motif = motif_file.read().strip()

        for seq in sequences:
            locations = [i for i in range(len(seq)) if seq[i:i+len(conserved_motif)] == conserved_motif]
            motif_locations.append(locations)
    else:
        common_subsequences = []
        min_length = min(len(seq) for seq in sequences)

        for length in range(min_length, 0, -1):
            for i in range(min_length - length + 1):
                subsequence = sequences[0][i:i+length]
                if all(subsequence in seq for seq in sequences[1:]):
                    common_subsequences.append(subsequence)

        motif_locations.append(common_subsequences)

    return motif_locations
