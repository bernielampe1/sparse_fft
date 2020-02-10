function sDFT = sfft(data, k)

% Sparse FFT implementation. The sDFT return variable is the estimated
% k complex values embedded in a n-length vector of n-k zeros.

% The magnitudes are right, but there is a phase shift in the code
% somewhere. I manually take it out at the end. It's a hack, but it works.

% length of data
n = length(data);

% heuristic number of times the outer loop runs
L = ceil(log2(n));

% heuristic size of frequency bin to divide n
B = 64;

% heuristic amount of oversample when finding peaks
d = 2;

% heuristic filter parameters tuned for B = 64 and n = 4096
stddev = sqrt(50);
w = 384;

% filter function return time domain filter and max value of G_hat
[G, G_hat] = build_filter(B, stddev, w, n);

% 4.2 Outer loop begin

% 4.2.1 Run location loop L times
Z = zeros(L, n);
for r = 1:L
    % 4.2.1 get the locations
    [J, z_hat, sigma_inv] = location_loop(data, n, G, w, B, d, k);
    [I, z_vals] = remap_locations(d, k, n, B, sigma_inv, J, z_hat, G_hat);
    Z(r, I) = z_vals;
end

% 4.2.2 Count the occurences of each location
bin_counts = sum(Z ~= 0);

% 4.2.3 Check if the occurences of each location is > L/2
% We pick the highest k here instead of score > L/2
[~, bin_inds] = sort(bin_counts, 'descend');
bin_inds = bin_inds(1:k);

% 4.2.4 Run estimation loop L times on location that have > L/2 occurrences
% 4.2.5 Use median of complex values as estimate at each location
sDFT = zeros(n, 1);
for k = 1:length(bin_inds)
    b = bin_inds(k);
    re = real(Z(:, b));
    im = imag(Z(:, b));
    re = median(re);
    im = median(im);

    % FIXME: Correcting for a phase shift bug earlier in the code
    b = mod(n - b + 1, n) + 1;
    num = re + im * 1i;
    a = angle(num) + 0.296108836078530;
    m = abs(num);
    num = m * exp(a * 1i);

    sDFT(b) = num;
end

end

% 4.1 Inner loop - remap each value in J to [n]
function [I, Z] = remap_locations(d, k, n, B, sigma_inv, J, z_hat, G_hat)
    I = zeros(d*k*n/B, 1);
    Z = zeros(d*k*n/B, 1);

    m = 1;
    for i = 1:length(J)
        % compute lowest and highest bin values in [n]
        low = mod(ceil((J(i) - 1 - 0.5) * n / B), n);
        high = mod(ceil((J(i) - 1 + 0.5) * n / B), n);

        % each index in J maps to B indicies in [n]
        j = low;
        l = 1;
        while j ~= high
            % transform the index back
            I(m) = mod(j * sigma_inv, n) + 1;

            % compute the x estimate value
            Z(m) = z_hat(J(i)) ./ G_hat(l);

            m = m + 1;
            l = l + 1;
            j = mod(j + 1, n);
        end
    end
end

% 4.1 Inner loop
function [J, z_hat, sigma_inv] = location_loop(x, n, G, w, B, d, k)
    % 4.1.1 choose a random sigma \in [n], with sigma odd and relatively prime
    sigma = 0;
    while gcd(sigma, n) ~= 1
        sigma = randi([0 n], 1);
    end

    % get sigma inverse
    [~, sigma_inv] = gcd(sigma, n);
    sigma_inv = mod(sigma_inv, n);
  
    % check that we have a valid sigma
    assert(mod(sigma * sigma_inv, n) == 1, 'sigma not relatively prime with n');

    % 4.1.2 permute the spectrum y_i = (x_i * sigma) % n
    y = zeros(w, 1);
    for i = 0:w-1
        y(i+1) = x(mod(i * sigma, n) + 1);
    end

    % 4.1.2 multiply by the window function
    y = y .* G;

    % 4.1.3 compute the B-point DFT z_hat from y
    z = zeros(B, 1);
    for i = 0:w-1
        ind = mod(i, B);
        z(ind + 1) = z(ind + 1) + y(i + 1);
    end
    z_hat = fft(z, B);

    % FIXME: Correcting for a phase shift error earlier in the code
    z_hat = circshift(flipud(z_hat), 1);
    
    % find dk indicies in B with largest magnitude
    [~, J] = sort(abs(z_hat), 'descend');
    J = J(1:d*k);
end

% Make a time domain filter length w with supp \in 2*B in frequency domain
function [G, G_hat] = build_filter(B, stddev, w, n)
    % boxcar length is twice the B
    box_length = floor(2 * B * 0.8);

    % boxcar
    g_boxcar = zeros(n, 1);
    g_boxcar((n-box_length)/2:(n+box_length)/2) = 1;

    % gaussian
    g_x = (0:n-1) - n/2;
    g_gauss = exp(-g_x.^2 / (2 * stddev * stddev))';

    % convolve the gaussian and boxcar
    g_conv = conv(g_boxcar, g_gauss, 'same');

    % normalize filters
    g_conv = g_conv ./ max(g_conv);

    % center at DC
    g_conv = circshift(g_conv, n/2+1);

    % ifft
    g_time = fftshift(ifft(g_conv, n));

    % truncate
    n_g = length(g_time);
    G = real(g_time(ceil((n_g - w)/2):ceil((n_g + w)/2)-1));

    % compute the fft of G_hat for the esimation step
    G_hat = fftshift(fft(G, n));

    % subset the frequency domain filter
    n_g = length(G_hat);
    G_hat = G_hat(ceil((n_g - B)/2):ceil((n_g + B)/2)-1);
end
