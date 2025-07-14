function K = gaussian_smoothing_matrix(N, kernel_size, sigma)
%GAUSSIAN_SMOOTHING_MATRIX Efficient NxN Gaussian smoothing matrix using toeplitz

    if mod(kernel_size, 2) == 0
        error('Kernel size must be odd');
    end

    % Generate centered 1D Gaussian kernel
    g = gaussian_kernel_1d(kernel_size, sigma);  % column vector

    % Embed kernel into a row vector of length N (zero-padded)
    pad = (N - kernel_size);
    padded = [zeros(1, floor(pad/2)), g', zeros(1, ceil(pad/2))];

    % Circularly shift to align kernel center at position 1
    center_index = find(g == max(g), 1);  % find peak (center)
    shift = center_index - 1;
    padded = circshift(padded, -shift);

    % Build Toeplitz matrix from padded kernel
    first_col = [padded(1), fliplr(padded(2:end))];
    first_row = padded;

    K = toeplitz(first_col, first_row);
end
