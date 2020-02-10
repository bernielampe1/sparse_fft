function ERR = l1_vs_k(maxk)

% size of data
n = 4096;

% allocate
ERR = zeros(maxk, 1);

% loop over k
for k = 1:maxk
    % create data set
    data_freq = zeros(n, 1);
    
    % create k magnitudes at different locations
    for j = 1:k
        data_freq(randi([1 4096])) = 1; %randi([1 100]);
    end
    
    % create time domain data
    data_time = ifft(data_freq) .* n;

    % run sparse fft
    dft = sfft(data_time, k);
    
    % compute l1 error
    e = sum(abs(dft - data_freq)) ./ k;

    ERR(k) = e;
end

end
