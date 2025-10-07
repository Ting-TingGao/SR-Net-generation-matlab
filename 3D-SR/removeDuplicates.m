function unique_matrix = removeDuplicates(matrix)
    % Sort each row in ascending order
    sorted_matrix = sort(matrix, 2);

    % Find the unique rows in the sorted matrix
    [~, unique_idx] = unique(sorted_matrix, 'rows', 'stable');

    % Extract the corresponding rows from the original matrix
    unique_matrix = matrix(unique_idx, :);
end