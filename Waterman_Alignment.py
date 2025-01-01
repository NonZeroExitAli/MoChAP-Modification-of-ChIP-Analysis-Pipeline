def read_fasta(file_path):
    headers = []
    sequences = []

    with open(file_path, 'r') as file:
        current_header = None
        current_sequence = ""

        for line in file:
            line = line.strip()
            # Use > or @ to generalize on FASTA files and FASTAQ files
            if line.startswith('>') or line.startswith("@"):
                if current_header is not None:
                    # Save the current header and sequence
                    headers.append(current_header)
                    sequences.append(current_sequence)
                    # Reset for the next sequence by separating by header
                    current_header = line[1:]
                    current_sequence = ""
                else:
                    # Write the first header
                    current_header = line[1:]
            else:
                current_sequence += line

        # Add the last sequence
        if current_header is not None:
            headers.append(current_header)
            sequences.append(current_sequence)

    return headers, sequences


def smith_waterman(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-2):
    n = len(seq1)
    m = len(seq2)

    score_matrix = [[0] * (m + 1) for _ in range(n + 1)]

    traceback_matrix = [[0] * (m + 1) for _ in range(n + 1)]

    max_score = 0
    max_i, max_j = 0, 0

    # Fill in the score matrix and traceback matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty

            # Determine the maximum score and update the traceback matrix
            current_score = max(0, match, delete, insert)
            score_matrix[i][j] = current_score

            if current_score > max_score:
                max_score = current_score
                max_i, max_j = i, j

            if current_score == match:
                traceback_matrix[i][j] = 0
            elif current_score == delete:
                traceback_matrix[i][j] = 1
            else:
                traceback_matrix[i][j] = 2

    # Traceback to find the aligned sequences
    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = max_i, max_j
    reference_start = j

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        if traceback_matrix[i][j] == 0:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 1:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1

    return aligned_seq1, aligned_seq2, reference_start



