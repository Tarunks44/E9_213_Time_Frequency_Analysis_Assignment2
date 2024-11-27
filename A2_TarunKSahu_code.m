%% Question 1: Generation of Discrete Time Sinusoid and Compute its Fourier Spectrum
% Parameters
% Need to change x axis of magnitude and phase from (-0.1 to 0.1)
close all;clear;clc
N = 512;  % Signal length
w0 = 20*pi/(2*N);  % Frequency in radians

% Time vector
n = -N:N;

% Generate the sinusoidal signal
s = cos(w0 * n);

% Compute the Discrete Fourier Transform (DFT) using custom function
S = my_fft(s);

% Frequency vector in radians per sample for plotting
omega = linspace(-pi, pi, length(S));

% Shift zero frequency component to the center using custom function
S_shifted = my_fftshift(S);

% Magnitude of the Fourier spectrum
S_mag = abs(S_shifted);

% Phase of the Fourier spectrum
S_phase = angle(S_shifted);

% Plot the time-domain signal
figure;
subplot(3, 1, 1);
stem(n, s, 'filled');
xlabel('n');
ylabel('s[n]');
title('Discrete-time Sinusoidal Signal');
grid on;

% Plot the magnitude of the Fourier spectrum with respect to omega
subplot(3, 1, 2);
stem(omega, S_mag);
xlabel('\omega (radians/sample)');
ylabel('|S(\omega)|');
title('Magnitude of Discrete Fourier Spectrum');
grid on;axis([-3 3 0 550])

% Plot the phase of the Fourier spectrum with respect to omega
subplot(3, 1, 3);
plot(omega, S_phase);
xlabel('\omega (radians/sample)');
ylabel('Phase (radians)');
title('Phase of Discrete Fourier Spectrum');
grid on;axis([-3 3 -4 4])

%% Question 2
% Generate a Gaussian signal with variance σ2 and compute its spectrum.
% Estimate σt and σω from the discrete Gaussian and its spectrum. 
% Compare σtσω with the theoretical lower limit.
% Parameters

% Parameters
N = 1024;  % Number of samples
sigma_t = 1;  % Time domain standard deviation
dt = 0.01;  % Time step

% Generate time vector
t = (-N/2:N/2-1) * dt;

% Generate Gaussian signal
x = exp(-t.^2 / (2*sigma_t^2)) / (sigma_t * sqrt(2*pi));

% Compute FFT and frequency vector
X = my_fftshift(my_fft(x)) * dt;
df = 1 / (N*dt);
f = (-N/2:N/2-1) * df;
omega = 2 * pi * f;

% Calculate time domain mean and standard deviation
t_mean = sum(t .* x) / sum(x);
sigma_t_calc = sqrt(sum((t - t_mean).^2 .* x) / sum(x));

% Calculate frequency domain mean and standard deviation
omega_mean = sum(omega .* abs(X).^2) / sum(abs(X).^2);
sigma_omega_calc = sqrt(sum((omega - omega_mean).^2 .* abs(X).^2) / sum(abs(X).^2));

% Calculate the product of standard deviations
time_freq_product = sigma_t_calc * sigma_omega_calc;

% Theoretical lower limit
lower_limit = 0.5;

% Create figure with 3 subplots
figure('Position', [100, 100, 800, 800]);

% Plot Gaussian signal
subplot(3,1,1);
plot(t, x);
title('Gaussian Signal in Time Domain');
xlabel('Time');
ylabel('Amplitude');
grid on

% Plot magnitude spectrum
subplot(3,1,2);
plot(omega, abs(X));
title('Magnitude Spectrum in ω Domain');
xlabel('Angular Frequency (ω)');
ylabel('|X(ω)|');
grid on

% Display results in text box
subplot(3,1,3);
axis off;
results_text = sprintf(['Time domain standard deviation (σ_t): %.4f\n' ...
                        'Frequency domain standard deviation (σ_ω): %.4f\n' ...
                        'Product σ_t * σ_ω: %.4f\n' ...
                        'Theoretical lower limit: %.4f\n'], ...
                        sigma_t_calc, sigma_omega_calc, time_freq_product, ...
                        lower_limit);
text(0.1, 0.5, results_text, 'FontSize', 12, 'VerticalAlignment', 'middle');
title('Numerical Results');


%% Question 3
% Generate Gabor functions and determine their spectra for varying ⟨t⟩, ⟨ω⟩ and σ2.


% Parameters
sigma2_values = [5, 20];  % Variance of the Gaussian window (reduce number of values)
t_mean_values = [-10, 0]; % Different time window means (reduce number of values)
omega_mean_values = [pi/4, pi/2];  % Different center frequencies (reduce number of values)

N = 1024;  % Number of points
t = linspace(-50, 50, N);  % Time vector

% Number of subplots needed
num_plots = length(sigma2_values) * length(t_mean_values) * length(omega_mean_values);

figure;

% Loop over different parameter values
plot_idx = 1;

for sigma2 = sigma2_values
    for t_mean = t_mean_values
        for omega_mean = omega_mean_values
            
            % Generate the Gabor function
            sigma = sqrt(sigma2);
            g = exp(-(t - t_mean).^2 / (2 * sigma2)) .* cos(omega_mean * (t - t_mean));
            
            % Compute the Fourier Transform
            G = my_fft(g);
            
            % Frequency vector in radians per sample
            omega = linspace(-pi, pi, N);
            
            % Shift zero frequency component to the center
            G_shifted = my_fftshift(G);
            
            % Magnitude of the Fourier spectrum
            G_mag = abs(G_shifted);
            
            % Plot the Gabor function in the time domain
            subplot(num_plots, 2, 2 * plot_idx - 1);
            plot(t, g);
            xlabel('Time (t)');
            ylabel('g(t)');
            title(sprintf('\\sigma^2 = %d, \\langle t \\rangle = %d, \\langle \\omega \\rangle = %.2f', sigma2, t_mean, omega_mean));
            grid on;
            
            % Plot the magnitude of the Fourier spectrum
            subplot(num_plots, 2, 2 * plot_idx);
            plot(omega, G_mag);
            xlabel('\omega (radians/sample)');
            ylabel('|G(\omega)|');
            title('Magnitude of Fourier Spectrum');
            grid on;
            
            plot_idx = plot_idx + 1;
        end
    end
end

%% Question No 4
% Generate a discrete-time Gaussian modulated chirp signal
% Parameters
% Parameters
N = 1024;  % Number of points
t = 0:N-1; % Time vector
sigma2 = 10000;  % Variance of the Gaussian envelope
fc_initial = 0.1;  % Initial frequency of the chirp
fc_final = 0.3;    % Final frequency of the chirp

% Generate the Gaussian modulated chirp signal
sigma = sqrt(sigma2);
chirp_signal = exp(-(t - N/2).^2 / (2 * sigma2)) .* cos(2 * pi * (fc_initial * t + (fc_final - fc_initial) * (t.^2) / (2 * N)));

% Compute the Fourier Transform
S = my_fft(chirp_signal);

% Frequency vector
f = linspace(-0.5, 0.5, N);

% Shift zero frequency component to the center
S_shifted = my_fftshift(S);

% Magnitude of the Fourier spectrum
S_mag = abs(S_shifted);

% Calculate the instantaneous frequency (f vs n)
instantaneous_frequency = fc_initial + (fc_final - fc_initial) * t / (N-1);

% Plot the f vs n (instantaneous frequency)
figure;
subplot(3, 1, 1);
plot(t, instantaneous_frequency, 'k-', 'LineWidth', 1.5);
xlabel('n →');
ylabel('frequency →');
title('f vs n Linear Plot');
grid on;
hold on;
plot([1024 1024], [0 0.3], 'k--'); % Line at N=1024 and f=0.3
hold off;

% Plot the chirp signal
subplot(3, 1, 2);
plot(t, chirp_signal);
xlabel('n →');
ylabel('Amplitude');
title('Gaussian Modulated Chirp Signal');
grid on;

% Plot the magnitude of the Fourier spectrum
subplot(3, 1, 3);
plot(f, S_mag);
xlabel('Normalized Frequency');
ylabel('Magnitude');
title('Magnitude Spectrum of the Chirp Signal');
grid on;


%% Question 5: Chirp Signal with Custom Function
% Parameters
N = 1024;  % Number of points
t = linspace(-10, 10, N); % Time vector

% Different choices of alpha, beta, and w0
alpha_values = [0.5, 1];
beta_values = [0.1, 0.5];
w0_values = [0, pi/4, pi/2];

% Create a figure for the plots
figure;

plot_idx = 1;

% Keep track of the number of plots
num_plots = 5;  % We will plot 5 different combinations of parameters

% Loop over different parameter combinations
for alpha = alpha_values
    for beta = beta_values
        for w0 = w0_values
            if plot_idx > num_plots
                break;
            end

            % Generate the function f(t) according to the given formula
            f_t = (alpha / pi)^(1/4) * exp(-alpha * t.^2 / 2 + 1i * beta * t.^2 / 2 + 1i * w0 * t);
            
            % Compute the Fourier Transform
            F = my_fft(f_t);
            
            % Frequency vector
            f = linspace(-0.5, 0.5, N);
            
            % Shift zero frequency component to the center
            F_shifted = my_fftshift(F);
            
            % Magnitude of the Fourier spectrum
            F_mag = abs(F_shifted);
            
            % Plot the time-domain function f(t)
            subplot(2, 5, plot_idx);
            plot(t, real(f_t));
            xlabel('t');
            ylabel('Re[f(t)]');
            title(sprintf('\\alpha = %.2f, \\beta = %.2f, w_0 = %.2f', alpha, beta, w0));
            grid on;
            
            % Plot the magnitude of the Fourier spectrum
            subplot(2, 5, plot_idx + 5);
            plot(f, F_mag);
            xlabel('Frequency');
            ylabel('|F(f)|');
            title('Magnitude Spectrum');
            grid on;
            
            plot_idx = plot_idx + 1;  % Move to the next subplot pair
        end
    end
end

% Adjust the layout for better readability
sgtitle('f(t) and its Spectrum for Different \alpha, \beta, and w_0'); % Super title
