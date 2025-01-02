import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from your_existing_file import SAM_Read, calculate_frequencies, calculate_motif_matrix, conserved_motif_file_exists, check_conserved_motif

def visualize_motif_matrix(motif_matrix, nucl):
    plt.figure(figsize=(8, 6))
    plt.imshow(motif_matrix['forward'], cmap='viridis', interpolation='nearest', aspect='auto')
    plt.title('Motif Matrix - Strand: FORWARD')
    plt.xlabel('nucl')
    plt.ylabel('Position')
    plt.xticks(np.arange(len(nucl)), list(nucl))
    plt.yticks(np.arange(len(motif_matrix['forward'])), np.arange(1, len(motif_matrix['forward']) + 1))
    plt.colorbar(label='Probability')
    plt.show()

def visualize_motif_locations(motif_locations):
    plt.figure(figsize=(10, 6))
    for idx, positions in enumerate(motif_locations, start=1):
        plt.scatter(positions, [idx] * len(positions), label=f'Sequence {idx}', marker='o', alpha=0.7)
    plt.title('Motif Locations in Sequences')
    plt.xlabel('Position in Sequence')
    plt.ylabel('Sequence Index')
    plt.yticks(np.arange(1, len(motif_locations) + 1))
    plt.legend()
    plt.grid(True)
    plt.show()

def visualize_venn_diagram(original_motif_locations, online_motif_locations):
    original_motif_set = set(tuple(loc) for loc in original_motif_locations)
    online_motif_set = set(tuple(loc) for _, loc in online_motif_locations)

    venn_labels = {'100': len(original_motif_set - online_motif_set),
                   '010': len(online_motif_set - original_motif_set),
                   '110': len(original_motif_set.intersection(online_motif_set))}

    plt.figure(figsize=(8, 8))
    venn_diagram = venn2(subsets=(len(original_motif_set), len(online_motif_set), len(original_motif_set.intersection(online_motif_set))),
                         set_labels=('Local Motif', 'Online Motif'))

    for text in venn_diagram.set_labels:
        text.set_fontsize(16)
    for text in venn_diagram.subset_labels:
        text.set_fontsize(14)

    plt.title('Venn Diagram - Motif Locations', fontsize=18)
    plt.show()

if __name__ == "__main__":
    SAM_file = "header.sam"
    original_sequences = SAM_Read(SAM_file)
    print("Original Sequences:")
    print(original_sequences)

    motif_file_path = "motif.txt"
    with open(motif_file_path, "r") as motif_file:
        motif = motif_file.read().strip()

    motif_locations = check_conserved_motif(original_sequences, motif)

    visualize_motif_locations(motif_locations)

    online_sequences_url = "https://your.online.source/sequences.txt"
    online_sequences = fetch_sequences_online(online_sequences_url)

    common_locations = compare_motif_locations_with_online_sequences(motif_locations, online_sequences)

    print("\nCommon Motif Locations:")
    for seq, loc in common_locations:
        print(f"Sequence: {seq}, Location: {loc}")

    motif_matrix = calculate_motif_matrix(original_sequences)
    print_motif_results(motif_matrix, original_sequences)

    visualize_motif_matrix(motif_matrix, set(''.join(original_sequences)))

    similarity_percentage = calculate_similarity_percentage(motif_locations, common_locations)
    print(f"\nSimilarity Percentage: {similarity_percentage:.2f}%")

    visualize_venn_diagram(motif_locations, common_locations)
