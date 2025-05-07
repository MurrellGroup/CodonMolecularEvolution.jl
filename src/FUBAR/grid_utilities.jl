"""
    generate_ambient_indices(N::Int)

Generate ambient format indices for an NxN grid. Indices are generated along diagonals
from bottom-right to top-left, with each diagonal numbered from bottom to top.
"""
function generate_ambient_indices(N::Int)
    indices = zeros(Int, N, N)
    current_index = 1

    # Loop over each diagonal, starting from bottom right
    for diag in (2N-2):-1:0
        # Calculate range for this diagonal
        col_start = max(0, diag - N + 1)
        col_end = min(diag, N - 1)

        # Store positions for this diagonal (bottom to top)
        positions = [(diag - col, col)
                    for col in col_end:-1:col_start]

        # Fill in indices
        for pos in positions
            indices[pos[1]+1, pos[2]+1] = current_index
            current_index += 1
        end
    end

    return indices
end

"""
    generate_fubar_indices(N::Int)

Generate FUBAR format indices for an NxN grid. Indices are generated column-wise
from bottom to top, starting from the leftmost column.
"""
function generate_fubar_indices(N::Int)
    indices = zeros(Int, N, N)
    current_index = 1

    # Fill column by column, bottom to top
    for col in 1:N
        for row in N:-1:1
            indices[row, col] = current_index
            current_index += 1
        end
    end

    return indices
end

function get_ambient_to_fubar_permutation(N::Int)
    ambient = generate_ambient_indices(N)
    fubar = generate_fubar_indices(N)

    # Create mapping from ambient indices to FUBAR indices
    n_elements = N * N
    perm = zeros(Int, n_elements)

    for i in 1:N
        for j in 1:N
            ambient_idx = ambient[i, j]
            fubar_idx = fubar[i, j]
            # We want to map FROM ambient TO fubar, so:
            perm[fubar_idx] = ambient_idx  # This line was wrong before!
        end
    end

    return perm
end

function get_fubar_to_ambient_permutation(N::Int)
    ambient = generate_ambient_indices(N)
    fubar = generate_fubar_indices(N)

    # Create mapping from FUBAR indices to ambient indices
    n_elements = N * N
    perm = zeros(Int, n_elements)

    for i in 1:N
        for j in 1:N
            ambient_idx = ambient[i, j]
            fubar_idx = fubar[i, j]
            # We want to map FROM fubar TO ambient, so:
            perm[ambient_idx] = fubar_idx  # This line was wrong before!
        end
    end

    return perm
end
