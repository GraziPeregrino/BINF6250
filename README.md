# Introduction
This project implements distance-based phylogenetic tree construction using the neighbor-joining (NJ) algorithm, with pairwise distances derived from Smith–Waterman local alignment scores. It extends our Project 05 alignment code into a full pipeline that takes a set of sequences and produces an unrooted tree in Newick format, along with a visualization.

The pipeline has three stages:
1. Distance calculation. Compute all pairwise Smith–Waterman alignment scores, convert similarity scores to distances, and assemble the complete distance matrix. A normalized similarity matrix can optionally be built as well.
2. Tree construction. Build the neighbor-joining matrix of rate-adjusted distances, then iteratively join the pair selected by the NJ criterion. After each join, update the distances with the NJ formula and compute branch lengths for the newly joined nodes.
3. Tree visualization. Convert the resulting tree structure to Newick format and generate an unrooted tree visualization.

# Pseudocode
```
FUNCTION READ_FASTA(filename):

    # INITIALIZATION
    sequences = EMPTY_DICTIONARY()
    current_id = NONE
    current_sequence = []

    # OPEN FASTA FILE
    OPEN filename FOR READING AS file

    # READ FILE LINE BY LINE
    FOR EACH line IN file:

        # Remove whitespace and newline characters
        line = STRIP(line)

        # Skip empty lines
        IF line == "":
            CONTINUE

        # CHECK FOR A NEW SEQUENCE HEADER
        IF line STARTS WITH ">":

            # Save the previous sequence if one exists
            IF current_id IS NOT NONE:
                sequences[current_id] = JOIN(current_sequence)

            # Extract sequence identifier
            current_id = FIRST_WORD(line WITHOUT FIRST CHARACTER)

            # Reset sequence list for new entry
            current_sequence = []

        ELSE:

            # Append sequence data
            APPEND line TO current_sequence

    # SAVE THE LAST SEQUENCE
    IF current_id IS NOT NONE:
        sequences[current_id] = JOIN(current_sequence)

    # CLOSE FILE
    CLOSE file

    RETURN sequences


FUNCTION SMITH_WATERMAN(sequence1, sequence2,
                        match_score = 2.0,
                        mismatch_penalty = -1.0,
                        gap_penalty = -1.0):

    n = LENGTH(sequence1)
    m = LENGTH(sequence2)

    # Handle empty sequences
    IF n == 0 OR m == 0:
        RETURN 0.0

    # Validate scoring parameters
    IF match_score <= 0 OR mismatch_penalty > 0 OR gap_penalty > 0:
        RAISE ERROR "Invalid scoring parameters"

    # Initialize matrices with zeros
    INITIALIZE score_matrix of size (n + 1) × (m + 1)
        filled with 0.0

    INITIALIZE traceback_matrix of size (n + 1) × (m + 1)
        filled with END (0)

    max_score = 0.0
    maximum_position = NONE

    # FORWARD SWEEP: Calculate scores and directions
    FOR i FROM 1 TO n:
        FOR j FROM 1 TO m:

            # Compare characters
            IF sequence1[i - 1] == sequence2[j - 1]:
                substitution = match_score
            ELSE:
                substitution = mismatch_penalty

            # Calculate possible scores
            diagonal = score_matrix[i - 1][j - 1] + substitution
            up = score_matrix[i - 1][j] + gap_penalty
            left = score_matrix[i][j - 1] + gap_penalty

            # Apply Smith-Waterman recurrence
            score = MAX(0, diagonal, up, left)

            score_matrix[i][j] = score

            # Record traceback direction
            IF score == 0:
                traceback_matrix[i][j] = END (0)

            ELSE IF score == diagonal:
                traceback_matrix[i][j] = DIAG (1)

            ELSE IF score == up:
                traceback_matrix[i][j] = UP (2)

            ELSE:
                traceback_matrix[i][j] = LEFT (3)

            # Track maximum score and its position
            IF score > max_score:
                max_score = score
                maximum_position = (i, j)

    # Handle sequences with no positive alignment
    IF max_score == 0:
        RETURN 0.0

    # BACKWARD SWEEP: Reconstruct and score alignment
    (aligned_seq1, aligned_seq2, alignment_score) =
        TRACEBACK_AND_SCORE(
            sequence1,
            sequence2,
            traceback_matrix,
            maximum_position,
            match_score,
            mismatch_penalty,
            gap_penalty
        )

    # Calculate self-alignment scores
    # Valid for positive match and nonpositive penalties
    score_uu = n * match_score
    score_vv = m * match_score

    # Normalize using assignment formula
    normalized_score = alignment_score / SQRT(score_uu * score_vv)

    RETURN normalized_score


FUNCTION TRACEBACK_AND_SCORE(seq1, seq2, traceback_matrix,
                             maximum_position,
                             match, mismatch, gap):

    # Initialize lists and alignment score
    aligned_seq1 = []
    aligned_seq2 = []
    alignment_score = 0.0

    (current_row, current_col) = maximum_position

    # BACKWARD SWEEP
    WHILE current_row > 0 AND current_col > 0:

        current_move = traceback_matrix[current_row][current_col]

        # Stop when reaching a zero-score cell
        IF current_move == END (0):
            BREAK

        # DIAGONAL: Match or mismatch
        ELSE IF current_move == DIAG (1):

            char1 = seq1[current_row - 1]
            char2 = seq2[current_col - 1]

            APPEND char1 TO aligned_seq1
            APPEND char2 TO aligned_seq2

            IF char1 == char2:
                alignment_score = alignment_score + match
            ELSE:
                alignment_score = alignment_score + mismatch

            current_row = current_row - 1
            current_col = current_col - 1

        # UP: Gap in sequence 2
        ELSE IF current_move == UP (2):

            APPEND seq1[current_row - 1] TO aligned_seq1
            APPEND "-" TO aligned_seq2

            alignment_score = alignment_score + gap

            current_row = current_row - 1

        # LEFT: Gap in sequence 1
        ELSE IF current_move == LEFT (3):

            APPEND "-" TO aligned_seq1
            APPEND seq2[current_col - 1] TO aligned_seq2

            alignment_score = alignment_score + gap

            current_col = current_col - 1

    # Reverse lists to recover alignment order
    aligned_seq1 = JOIN(REVERSE(aligned_seq1))
    aligned_seq2 = JOIN(REVERSE(aligned_seq2))

    RETURN aligned_seq1, aligned_seq2, alignment_score

FUNCTION BUILD_DISTANCE_MATRIX(sequences):

    # INITIALIZATION
    ids = LIST of keys from sequences
    n = LENGTH(ids)

    # Initialize n × n distance matrix with zeros
    Initialize distance_matrix of size n × n
    filled with 0.0

    # CALCULATE PAIRWISE DISTANCES
    FOR i FROM 0 TO n - 1:

        # Get first sequence
        seq1 = sequences[ids[i]]

        FOR j FROM i + 1 TO n:

            # Get second sequence
            seq2 = sequences[ids[j]]

            # Calculate normalized Smith-Waterman similarity
            similarity = SMITH_WATERMAN(seq1, seq2)

            # Convert similarity into distance
            distance = 1.0 - similarity

            # Update symmetric distance matrix
            distance_matrix[i][j] = distance
            distance_matrix[j][i] = distance

    # RETURN RESULTS
    RETURN distance_matrix, ids



FUNCTION NEIGHBOR_JOINING(distance_matrix, labels):

    # INITIALIZATION
    active_nodes = COPY(labels)
    n = LENGTH(active_nodes)

    IF n < 2:
        RAISE ERROR "At least two sequences are required"

    # Convert numerical matrix to node-based distances
    D = EMPTY_DICTIONARY()

    FOR i FROM 0 TO n - 1:
        D[labels[i]] = EMPTY_DICTIONARY()

        FOR j FROM 0 TO n - 1:
            D[labels[i]][labels[j]] = distance_matrix[i][j]

    # Store Newick string for each node
    newick = EMPTY_DICTIONARY()

    FOR EACH node IN active_nodes:
        newick[node] = node

    node_counter = 0

    # ITERATION: Continue until two nodes remain
    WHILE n > 2:

        # STEP 1: Calculate row sums
        row_sums = {}

        FOR EACH i IN active_nodes:
            row_sums[i] = 0

            FOR EACH j IN active_nodes:
                row_sums[i] = row_sums[i] + D[i][j]

        # STEP 2: Calculate Q-values and find minimum pair
        min_q = INFINITY
        min_pair = NONE

        FOR EACH unique pair (i, j) IN active_nodes:

            q = (n - 2) * D[i][j] -
                row_sums[i] - row_sums[j]

            IF q < min_q:
                min_q = q
                min_pair = (i, j)

        # STEP 3: Select pair
        (i, j) = min_pair

        # STEP 4: Calculate branch lengths
        branch_i = 0.5 * (
            D[i][j] +
            (row_sums[i] - row_sums[j]) / (n - 2)
        )

        branch_j = D[i][j] - branch_i

        # STEP 5: Create a new internal node
        k = GENERATE_UNIQUE_NODE_ID(node_counter)
        node_counter = node_counter + 1

        # Combine the two subtrees into a Newick string
        newick[k] = "(" +
            newick[i] + ":" + FORMAT(branch_i) + "," +
            newick[j] + ":" + FORMAT(branch_j) + ")"

        # STEP 6: Calculate distances to remaining nodes
        remaining_nodes = active_nodes EXCLUDING i AND j

        D[k] = EMPTY_DICTIONARY()

        FOR EACH x IN remaining_nodes:

            distance_kx = 0.5 * (
                D[i][x] + D[j][x] - D[i][j]
            )
            D[k][x] = distance_kx
            D[x][k] = distance_kx

        D[k][k] = 0

        # STEP 7: Remove joined nodes from active matrix
        REMOVE D[i] AND D[j]

        FOR EACH x IN remaining_nodes:
            REMOVE D[x][i] AND D[x][j]

        # Remove old Newick entries
        REMOVE newick[i] AND newick[j]

        # STEP 8: Update active nodes
        active_nodes = remaining_nodes
        APPEND k TO active_nodes

        n = LENGTH(active_nodes)

    # TERMINATION: Two nodes remain
    (i, j) = active_nodes

    final_distance = D[i][j]

    # Split final edge to represent the unrooted tree
    final_branch = final_distance / 2

    final_tree = "(" +
        newick[i] + ":" + FORMAT(final_branch) + "," +
        newick[j] + ":" + FORMAT(final_branch) + ");"

    RETURN final_tree

```

# Testing our code
We conducted a series of tests on our code throughout the process, to make sure things were working as we expected. 
To test `read_fasta()`, we conducted the following:
```
# Test read_fasta()
sequences = read_fasta("data/lafayette_SARS_RT.fasta")

print("Number of sequences:", len(sequences))

# Display the first 3 sequences
for seq_id, sequence in list(sequences.items())[:3]:
    print(f"\nSequence ID: {seq_id}")
    print(f"Length: {len(sequence)}")
    print(f"First 50 characters: {sequence[:50]}")
```
Result:
```
Number of sequences: 20

Sequence ID: gi|24209940|gb|AAN41462.1|
Length: 229
First 50 characters: PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKI

Sequence ID: gi|24209942|gb|AAN41463.1|
Length: 229
First 50 characters: PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKI

Sequence ID: gi|24209944|gb|AAN41464.1|
Length: 229
First 50 characters: PISPIETVPVKLKPGMDGPKVKQWPLTEEKIKALVEICTEMEKEGKISKI
```

To test our `smith_waterman()` implementation, we conducted the following:
```
# Test 1: Identical sequences
print(smith_waterman("ACGT", "ACGT"))
# Expected: 1.0

# Test 2: Partial similarity
print(smith_waterman("ACGT", "ACGA"))
# Expected: 0.75

# Test 3: No matching characters
print(smith_waterman("AAAA", "TTTT"))
# Expected: 0.0

# Test 4: Empty sequence
print(smith_waterman("", "ACGT"))
# Expected: 0.0

# Test 5: Different lengths
print(smith_waterman("ACGT", "TTACGTAA"))
# Expected: approximately 0.71
```
Result:
```
1.0
0.75
0.0
0.0
0.71
```

To test `build_distance_matrix()`, we conducted the follwoing:
```
sequences = read_fasta("data/lafayette_SARS_RT.fasta")

distance_matrix, ids = build_distance_matrix(sequences)

print("Matrix shape:", distance_matrix.shape)
print("Sequence count:", len(ids))

print("\nFirst three rows and columns:")
print(np.round(distance_matrix[:3, :3], 4))
```
Result:
```
Matrix shape: (20, 20)
Sequence count: 20

First three rows and columns:
[[0.   0.01 0.05]
 [0.01 0.   0.04]
 [0.05 0.04 0.  ]]
```

# Successes
Description of the team's learning points

# Struggles
**Concept:** We initially struggled with the main idea behind this program, and how we needed to go about implmenting it. We spent a couple of ours of a Teams call trying to figure out if/how we needed to change our `smith_waterman()` implementation to get the results we needed. Eventually, we decided on slightly changing the implementations from Graziano's group from last project, to include the overall similarity score of the two sequences.

# Personal Reflections
## Group Leader - Graziano 
Group leader's reflection on the project

## Member - Margareth  
Other members' reflections on the project

# Generative AI Appendix
1. Description of which generative AI was used and its version.

NotebookLM 2.0 and Claude Opus 5 

2. The entire prompt that was used to generate the content.

Can you help me to fix my pseudocode to a the github format?
I need explanation in a different way of how the distance metrics for comparing sequences using Smith-Waterman avoiding the use of UPGMA
Can you help to understand the Neighbor - Joining algorithm and the use of the algorithm in bioinformatics

3. An explanation of how it was used (e.g., to generate ideas).
I used it as a tool to learn the concepts described on the articles and materials, asking for further clarification on the material and extra examples.

4. A justification for why generative AI was used.

Sometimes just reading a concept and formula is not enough for me to have a grasp on the material, so I used it as an additional source to lean about the Smith-Waterman algorithm, the use of traceback and scoring and the application of the Neighbor- Joining Algorithm
