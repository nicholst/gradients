function output_matrix = keep_top(input_matrix, r)
    % Input: 
    %   input_matrix - input matrix 
    %   r - fraction of top elements to keep (between 0 and 1)
    
    % Validate input
    if r < 0 || r > 1
        error('Fraction must be between 0 and 1');
    end
    
    % Create output matrix as a copy of input matrix
    output_matrix = input_matrix;
    
    % Iterate through each row
    for i = 1:size(input_matrix, 1)
        % Get current row
        row = input_matrix(i, :);
        
        % Calculate number of elements to keep
        num_keep = max(1, round(length(row) * r));
        
        % Sort row in descending order
        [sorted_row, indices] = sort(row, 'descend');
        
        % Zero out elements not in top fraction
        zero_mask = false(1, length(row));
        zero_mask(indices(1:num_keep)) = true;
        
        % Apply mask to output matrix row
        output_matrix(i, ~zero_mask) = 0;
    end
end