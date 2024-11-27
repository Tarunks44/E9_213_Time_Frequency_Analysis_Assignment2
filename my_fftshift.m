%% Custom FFTShift function
function X_shifted = my_fftshift(X)
    N = length(X);  % Length of the input spectrum
    X_shifted = zeros(1, N);  % Initialize the shifted result

    % Split the spectrum in half and swap
    if mod(N, 2) == 0
        % For even length signals
        mid = N / 2;
        X_shifted(1:mid) = X(mid+1:N);
        X_shifted(mid+1:N) = X(1:mid);
    else
        % For odd length signals
        mid = (N - 1) / 2;
        X_shifted(1:mid+1) = X(mid+1:N);
        X_shifted(mid+2:N) = X(1:mid);
    end
end