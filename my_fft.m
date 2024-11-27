%% Custom FFT function
function X = my_fft(x)
    N = length(x);  % Length of the input signal
    X = zeros(1, N);  % Initialize the result

    % Create the DFT matrix
    for k = 0:N-1
        for n = 0:N-1
            X(k+1) = X(k+1) + x(n+1) * exp(-1j * 2 * pi * k * n / N);
        end
    end
end
