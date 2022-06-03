%% Assignment 2, Part 2 (Choosing a landing site)
%  Do not change before line 32
%  You will need to have generated A2P2Data.mat from 
%  GenerateAssignment2Data.m before working with this file.

%  Clearing and preparing the workspace
clear; clc; close all;

%  Load assignment data from A2P3Data.mat.
load('A2P2Data.mat');  

%=================================================================%
%
% Refer to the assignment sheet for details.
% Names of the variables are important,
% e.g. 'a1' is considered a different variable to 'A1'.
% Make sure variables have been declared as they appear in the brief.
%
% t - time domain vector
% f - frequency vector
% SIG - Fourier Transform
% T - selected value
% Noisesig - estimated noise signal
% a0,an,bn - Trigonometric Fourier Series Coefficients
% OR
% c0,cn - Complex Fourier Series Coefficients
% Noisesig_fs - approximation of noise
% im1 - image 1
% im2 - image 2
% 
%====Enter your code below this line================================

%% 2.1 - View Noisy Image
w = 640;                                            % Width  [px]
h = 480;                                            % Height [px]
pixels = w*h;                                       % Total pixels
figure('Name', "Noisy Image")                   
imshow(reshape(sig(1,:), h, w))

%% 2.2 - Time and Frequency Reference Vectors
Fs = 1000;                                          % Pixel Sampling Rate [Hz]
T_image = pixels/Fs;                                % Time to recieve an image
t = linspace(0, T_image, length(sig)+1);            % Time vector for full image
f = linspace(-Fs/2, Fs/2, length(sig)+1);           % Frequency Vector
f(end) = []; t(end)=[];                             % Remove last elements

%% 2.3 - Visualise Time and Frequency Domain of First Image
figure('Position', [400, 400, 800, 400]); hold on;
subplot(1, 2, 1)
plot(t, sig(1, :));
xlim([0 3])
xlabel("Time [s]"), ylabel("Amplitude"), title("Time Domain (First 3s)")

SIG = fft(sig, length(sig), 2);             % Fourier Transform of each image
subplot(1, 2, 2)
plot(f, abs(fftshift(SIG(1, :)/Fs)))
xlabel("Frequency [Hz]"), ylabel("Magnitude"), title("Frequency Domain")

%% 2.4 - Estimate periodic noise
% Visually, noise has period of 1.477s
T = candidateT(1);                                                      % period of noise
periodInSamples = T * Fs;                                               % samples in noise period
Noisesig = estimateNoise(sig(1, :), periodInSamples);                   % estimated noise function
NoisesigFull = repmat(Noisesig, [1 ceil(length(t)/length(Noisesig))]);  % repeat noise function to fill time domain
NoisesigFull = NoisesigFull(1:length(t));                               % truncate noise function to match time domain
% Plot and compare the noise signal to the recieved signal
figure('Name', "Periodic Noise Comparison"),
subplot(2, 1, 1)
plot(t, sig(1, :));
xlim([0 3])
xlabel("Time [s]"), ylabel("Amplitude"), title("Time Domain (First 3s)")
subplot(2, 1, 2)
plot(t, NoisesigFull)
xlim([0 3])
xlabel("Time [s]"), ylabel("Amplitude"), title("Periodic Noise Signal (First 3s)")

%% 2.5 - Model periodic noise
% Complex Fourier Coefficients
t1= t(t<T); t1(end) = [];                           % One Period Time vector (0-1.477s)
ts = t(2)-t(1);                                     % Sampling interval [s]
f0 = 1/T;                                           % Fundamental Frequency
N = 6;                                              % Number of harmonics
n = (-N:N).';                                       % Vector of harmonics
cn = f0 * Noisesig * exp(-1j*2*pi*f0*n*t1).' * ts;  % Fourier coefficients
c0 = cn(N+1);                                       % DC coefficient

%% 2.6 - Remove DC bias
% Bias - make the real coefficient c0 = 0;
c0 = 0;                                             % Direct DC coefficient
cn(N+1) = c0;                                       % Corresponding cn coefficient

%% 2.7 - Generate Approximation
Noisesig_fs = cn * exp(1j*2*pi*f0*n*t);             % Fourier Series for full t

%% 2.8 - Compare Approximation to Noise for one period
figure('Name', "Fourier Approximation of Noisesig"), hold on;
plot(t1, Noisesig)
plot(t1, Noisesig_fs(1:length(t1)))
xlim([0 T])
title("Comparison of Fourier Series and Noisesig (1 Period)"), xlabel("Time [s]"), ylabel("Amplitude")
legend("Noisesig", "Noisesig\_fs")

% optionally Plot alternate harmonics
Noisesig_9 = complexFS(Noisesig, t1, 9, t);
Noisesig_10 = complexFS(Noisesig, t1, 10, t);
figure('Name', "Alternate Harmonics"), hold on;
subplot(1, 3, 1), hold on;
plot(t1, Noisesig)
plot(t1, Noisesig_fs(1:length(t1)))
xlim([0 T])
title("Fourier approximation (6 harmonics)"), xlabel("Time [s]"), ylabel("Amplitude")
legend("Noisesig", "Noisesig\_fs")
subplot(1, 3, 2), hold on;
plot(t1, Noisesig)
plot(t1, Noisesig_9(1:length(t1)))
xlim([0 T])
title("Fourier approximation (8 harmonics)"),xlabel("Time [s]"), ylabel("Amplitude")
legend("Noisesig", "Noisesig\_fs")
subplot(1, 3, 3), hold on;
plot(t1, Noisesig)
plot(t1, Noisesig_10(1:length(t1)))
xlim([0 T])
title("Fourier approximation (10 harmonics)"),xlabel("Time [s]"), ylabel("Amplitude")
legend("Noisesig", "Noisesig\_fs")
% Looks like some of the more detailed noise at the peaks and valleys are
% not modelled by only 6 harmonics. Might be worth using 10 to get these, as 10 gets the little bits heaps good

%% 2.9 - Denoise the image
Noisesig_fs = Noisesig_10;                          % Use 10 harmonics
im1(1,:) = sig(1,:) - Noisesig_fs;                  % Subtract noise from signal
 
figure('Name', "Denoised Image"), hold on;
subplot(1, 2, 1)
imshow(reshape(im1(1,:), h, w))                     % Show denoised image
title("Partially Denoised Image")
subplot(1, 2, 2)
plot(f, abs(fftshift(fft(im1(1, :)/Fs))))           % Plot frequency domain
title("Magnitude Spectrum of Partially Denoised Image"), xlabel("Frequency [Hz]"), ylabel("Magnitude")

%% 2.10 - Remove bandlimited noise

% % evaluate the fourier transform of im1
IM1(1,:) = fft(im1(1,:));                           % Fourier transform of denoised hiddenSignal
% figure('Name', "Magnitude Spectrum of Partially Denoised Image")
% plot(f, abs(fftshift(IM1(1, :)/Fs)))                % Plot frequency domain in band of interest
% xlim([215 250])
% xlabel("Freuqency [Hz]"), ylabel("Magnitude"), title("Magnitude spectrum of bandlimited noise")
% 
% % define Band-stop filter
B_low = 223;                                        % Lower bandwidth frequency
B_high = 240;                                       % Upper bandwidth frequency
filter = ones(size(f));                             % Initiallise with ones
filter((abs(f) > B_low) & (abs(f) < B_high)) = 0;   % Assign 0 to elements in band
filtered = fftshift(IM1(1,:)) .* filter;            % Filter shifted signal in frequency domain
im2 = zeros(size(sig));                             % Initialise for performance
im2(1,:) = ifft(ifftshift(filtered));               % Filtered image signal in time domain

figure('Name', "Filtered Image")
imshow(reshape(im2(1,:), h, w))
%% 2.11 - Repeat denoising
im1(2, :) = sig(2, :) - Noisesig_fs;                % Denoise image 2
im1(3, :) = sig(3, :) - Noisesig_fs;                % Denoise image 3
im1(4, :) = sig(4, :) - Noisesig_fs;                % Denoise image 4
IM1 = fft(im1, length(im1), 2);                     % Row-wise Fourier transform
% Determine bandwidth of bandlimited noise
figure('Name', "Magnitude Spectrum of Other Images with Bandlimited Noise"), hold on;
for k = 2:size(sig, 1)
    subplot(1, size(sig, 1)-1, k-1)
    plot(f, abs(fftshift(IM1(k, :)/Fs)))
    title(sprintf("Magnitude Spectrum of Image %i", k))
    xlabel("Frequency [Hz]"), ylabel("Magnitude")
end
% Visually get bands of noise
bands = {[223, 240], [232, 246], [234, 251], [253, 268]};
% Remove bandlimited noise from each image
filter = ones(size(f));                             % Initialise for performance
im2 = zeros(size(im1));                             % Initialise for performance
for k = 1:length(bands)
    filter(1:end) = 1;                              % Reallocate ones
    filter(abs(f) > bands{k}(1) & abs(f) < bands{k}(2)) = 0;
    filteredIm = fftshift(IM1(k, :)) .* filter;     % Apply the filter in frequency domain
    im2(k, :) = ifft(ifftshift(filteredIm));        % Store the final image in time domain
end
% Show all images
figure('Name', "All filtered Images"), hold on;
for k = 1:size(im2, 1)
    subplot(2, 2, k)
    imshow(reshape(im2(k, :), h, w))
    title(sprintf("Potential Landing Site %i", k))
end

%% 2.12


%% Helper functions
function [sig_fs] = complexFS(sig, t1, N, t)
    ts = t1(2)-t1(1);                               % Sampling interval [s]
    T = t1(end)-t1(1)+ts;                           % Period
    f0 = 1/T;                                       % Fundamental Frequency
    n = (-N:N).';                                   % Vector of harmonics
    cn = f0 * sig * exp(-1j*2*pi*f0*n*t1).' * ts;   % Fourier coefficients
    cn(N+1) = 0;                                    % Remove DC
    sig_fs = cn * exp(1j*2*pi*f0*n*t);              % Fourier Series for full t
end
