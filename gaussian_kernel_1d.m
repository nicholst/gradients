function k = gaussian_kernel_1d(size, sigma)
    if mod(size, 2) == 0
        error('Kernel size must be odd');
    end
    x = -floor(size/2):floor(size/2);
    k = exp(-x.^2 / (2 * sigma^2));
    k = k / sum(k);
end
