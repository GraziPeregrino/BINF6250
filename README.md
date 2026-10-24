# Introduction
Description of the project

# Pseudocode
```
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

        FOR j FROM i + 1 TO n - 1:

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


FUNCTION NEIGHBOR_JOINING(D, sequences):

    # INITIALIZATION
    active_nodes = LIST(sequences)
    n = LENGTH(active_nodes)

    # Initialize tree with one leaf per sequence
    T = EMPTY_TREE()

    FOR EACH node IN active_nodes:
        ADD_NODE(T, node)

    # ITERATION: Continue until two nodes remain
    WHILE n > 2:

        # STEP 1: Calculate row sums
        row_sums = {}

        FOR EACH i IN active_nodes:
            row_sums[i] = 0

            FOR EACH j IN active_nodes:
                row_sums[i] = row_sums[i] + D[i][j]

        # STEP 2: Calculate Q-matrix and find minimum pair
        min_q = INFINITY
        min_pair = NONE

        FOR EACH unique pair (i, j) IN active_nodes:

            q = (n - 2) * D[i][j] -
                row_sums[i] - row_sums[j]

            IF q < min_q:
                min_q = q
                min_pair = (i, j)

        # STEP 3: Select pair with minimum Q-value
        (i, j) = min_pair

        # STEP 4: Calculate branch lengths
        branch_i = 0.5 * (
            D[i][j] +
            (row_sums[i] - row_sums[j]) / (n - 2)
        )

        branch_j = D[i][j] - branch_i

        # STEP 5: Create new internal node
        k = GENERATE_UNIQUE_NODE()

        ADD_NODE(T, k)
        ADD_EDGE(T, k, i, branch_i)
        ADD_EDGE(T, k, j, branch_j)

        # STEP 6: Calculate distances to remaining nodes
        remaining_nodes = active_nodes EXCLUDING i AND j

        FOR EACH x IN remaining_nodes:

            distance_kx = 0.5 * (
                D[i][x] + D[j][x] - D[i][j]
            )

            D[k][x] = distance_kx
            D[x][k] = distance_kx

        D[k][k] = 0

        # STEP 7: Remove joined nodes from active matrix
        REMOVE rows and columns i, j FROM D

        # STEP 8: Update active nodes
        active_nodes = remaining_nodes
        APPEND k TO active_nodes

        n = LENGTH(active_nodes)

    # TERMINATION: Connect remaining two nodes
    (i, j) = active_nodes

    ADD_EDGE(T, i, j, D[i][j])

    RETURN T


```

# Successes
Description of the team's learning points

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
Other members' reflections on the project

# Generative AI Appendix
As per the syllabus
