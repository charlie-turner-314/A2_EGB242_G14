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

%% 2.1
w = 640;
h = 480;
pixels = w*h;
figure(1)
imshow(reshape(sig(1,:), h, w))

%% 2.2
Fs = 1000;                                  % Pixel Sampling Rate [Hz]
T = pixels/Fs;                              % Time to recieve an image
t = linspace(0, T, length(sig)+1);          % Time vector for full image
f = linspace(-Fs/2, Fs/2, length(sig)+1);   % Frequency Vector
f(end) = []; t(end)=[];

%% 2.3
figure(2); hold on;
subplot(3, 1, 1)
plot(t, sig(1, :));
xlim([0 3])
xlabel("Time [s]"), ylabel("Amplitude"), title("Time Domain (First 3s)")

SIG = fft(sig, length(sig), 2);
subplot(3, 1, 2)
plot(f, abs(fftshift(SIG(1, :)/Fs)))
xlabel("Frequency [Hz]"), ylabel("Magnitude"), title("Frequency Domain")

% Looks like the noise is about 1.48 secs, so use candidateT(1)

%% 2.4
noisePeriod = candidateT(1);
noisePeriodSamples = noisePeriod * Fs;
Noisesig = estimateNoise(sig(1, :), noisePeriodSamples);
NoisesigFull = repmat(Noisesig, [1 ceil(length(t)/length(Noisesig))]);
NoisesigFull = NoisesigFull(1:length(t));
subplot(3, 1, 3)
plot(t, NoisesigFull)
xlim([0 3])
xlabel("Time [s]"), ylabel("Amplitude"), title("Periodic Noise Signal (First 3 s)")

%% 2.5
% Complex Fourier Coefficients
tNoise = linspace(0, noisePeriod, length(Noisesig)+1); tNoise(end) = [];
ts = tNoise(2)-tNoise(1); 
f0 = 1/noisePeriod;
N = 15;
n = (-N:N).';
cn = f0 * Noisesig * exp(-1j*2*pi*f0*n*tNoise).' * ts;
c0 = cn(N+1);
NoisesigApprox = cn * exp(1j*2*pi*n*f0*tNoise);




%% 2.6
% Bias - make the real coefficient c0 = 0;
cn(N+1) = 0;

%% 2.7
Noisesig_fs = cn * exp(1j*2*pi*f0*n*tNoise);

%% 2.8
figure(3), hold on;
plot(tNoise, Noisesig)
plot(tNoise, NoisesigApprox)
plot(tNoise, Noisesig_fs)
xlim([0 noisePeriod])
legend("Noisesig", "NoisesigApprox", "Noisesig\_fs")
% Looks like some of the more detailed noise at the peaks and valleys are
% not modelled by only 6 harmonics. Might be worth using 10 to get these (9
% has better overall shape but 10 gets the little bits heaps good)

%% 2.9

t_tot = length(t)*ts / noisePeriod;                 % no. of periods in total
Noisesig_fs = repmat(Noisesig_fs, [1, ceil(t_tot)]);    % fourier approx of periodic noise signal repeated over total time (rounded up)
Noisesig_fs = Noisesig_fs(1:length(t));                % fourier approx repeated over total time (exact)

im1(1,:) = sig(1,:) - Noisesig_fs;                      % de-noised signal
 
figure(4)
imshow(reshape(im1(1,:), h, w))

%% 2.10

% evaluate the fourier transform of hdm
IM1(1,:) = fft(im1(1,:));                                                     %fourier transform of demodulated hiddenSignal

% define bandwidth
B_upper = 223;                                                              % upper noise bandwidth
B_lower = 240;                                                              % lower noise bandwidth
B_upper_i = B_upper*(length(sig(1,:))*ts);                                  %upper index value 
B_lower_i = B_lower*(length(sig(1,:))*ts);                                  %lower index value
middle_f = length(sig(1,:))/2;

% create vector to filter signal
filter = ones(1, length(f));                                                % ones vector
filter(middle_f - B_upper_i : middle_f - B_lower_i) = 0;                    % filter
filter(middle_f + B_lower_i : middle_f + B_upper_i) = 0;                    % filter

% filter HDM with filter vector k_filter
filtered = IM1(1,:).*filter;                                        %filtered demodulated signal in frequency domain

% compute inverse fourier transform to change to time domain
im2(1,:) = ifft(filtered);     %%filtered demodulated signal in frequency domain

figure(5)
imshow(reshape(im2(1,:), h, w))
%% 2.11
%% 2.12
