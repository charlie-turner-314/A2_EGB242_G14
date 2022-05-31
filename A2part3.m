%% Assignment 2, Part 3 (Impulse room signal analysis)
%  Do not change before line 35
%  You will need to have generated A2P3Data.mat from 
%  GenerateAssignment2Data.m before working with this file.

%  Clearing and preparing the workspace
clear; clc; close all;

%  Load assignment data from A2P1Data.mat.
load('A2P3Data.mat');  

%=================================================================%
%
% Refer to the assignment sheet for details.
% Names of the variables are important,
% e.g. 'a1' is considered a different variable to 'A1'.
% Make sure variables have been declared as they appear in the brief.
%
% Ts - Sampling period
% t - Time domain vector
% MUXSIG - Frequency domain representation of mux
% f - Frequency vector
% fshift - Frequency shifts
% Mag - Magnitude
% Phase - Phase
% xdm - All de-shifted signals in the time domain
% XDM - Frequency domain representation of xdm 
% freqResponse - Frequency response of systems 
% impresp - Impulse response of chosen system 
% MSG - Frequency domain representation of Filtered signals 
% imp - Recovered impulse responses resampled to 48kHz 
% text - Recovered text streams resampled to 48kHz 
% t_recov - New time vector for recovered and resampled signals 
% f_recov - New frequency vector for recovered and resampled signals 
%
%====Enter your code below this line================================

%% 3.1
Ts = 1/fs;                                          % Sampling Period [s]
t = 0:Ts:0.5; t(end) = [];       % Time vector
MUXSIG = fft(muxSignal);                            % Fourier Transform
f = linspace(-fs/2, fs/2, length(t)+1); f(end)=[];  % Freuquency Vector

% Plot signal against time
figure('Name',"Recieved Signal"), hold on;
subplot(1, 2, 1)
plot(t, muxSignal);
xlabel("Time [s]"), ylabel("Amplitude"), title("Time Domain of muxSignal")


% Plot f
subplot(1, 2, 2)
plot(f, abs(fftshift(MUXSIG/fs)));
xlabel("Frequency [Hz]"), ylabel("Magnitude"), title("Magnitude Spectrum of muxSignal")

%% 3.2

fshift = [32e3 80e3 128e3 176e3 224e3 272e3];       % Frequency shifts
fshift_indices = zeros(size(fshift));
MUXSIG_shift = fftshift(MUXSIG);
for k = 1:length(fshift)
    fshift_indices(k) = find(f == fshift(k));
end
Mag = abs(MUXSIG_shift(fshift_indices))/fs;         % magnitude
Phase = angle(MUXSIG_shift(fshift_indices));        % [rads]    

%% 3.3
xdm = FDMDemux(muxSignal, t, Mag, fshift, Phase);
figure, hold on;
rows = size(xdm, 1);
for k = 1:rows
    subplot(rows, 2, 2*k-1)
    plot(t, xdm(k, :));
end
%% 3.4
XDM = fft(xdm, length(xdm), 2);
for k = 1:rows
    subplot(rows, 2, 2*k)
    plot(f, abs(fftshift(XDM(k, :)/fs)));
    title(sprintf("Magnitude Spectrum of Stream %i", k))
end
%% 3.5
% bode plot shows low pass is number 1. Make some stuff ap about the
% denominator lol. Also pooles have to be neg apparently
factorTF(sys(1)) % shows its  a time multiplexed exponential -> (1/(s+a)^2)
figure, title("System 1 Freq Response")
h = bodeplot(sys(4));
setoptions(h, 'FreqUnits', 'Hz')
%% 3.6

%% 3.7
%% 3.8  don't know anything after this
impulse = [1 zeros([1 length(XDM)-1])];
impulse_response = lsim(sys(1), impulse, t);
FrequencyResponse = fft(impulse_response);
MSG = zeros(size(XDM));
for k = 1:size(XDM, 1)
    MSG(k, :) = FrequencyResponse.' .* XDM(k, :);
end
msg = ifft(MSG, length(MSG), 2);
for k = 1:size(msg,1)
    msg(k, :) = msg(k, :) - mean(msg(k,:)); % remove dc bias (assuming mean should be 0);
    msg(k, :) = msg(k, :) / (max(abs(msg(k, :))));
end

figure, hold on
for k = 1:size(MSG, 1)
    subplot(size(MSG, 1), 2, 2*k-1)
    plot(t, msg(k, :))
    title(sprintf("Time Domain of Stream %i", k)), xlabel("Time [s]"), ylabel("Amplitude")
    ylim([-1, 1])
    subplot(size(MSG, 1), 2, 2*k)
    plot(f, real(abs(fftshift(MSG(k, :)/fs))))
    title(sprintf("Magnitude Spectrum of Stream %i", k)), xlabel("Frequency [Hz]"), ylabel("Magnitude")
    xlim([-10e3 10e3])
end
%% 3.9
fs_recov = 16e3;
impulses_recov = resample(msg((1:2:end), :).', fs_recov, fs).';
%% 3.10
texts_recov = msg((2:2:end), :);
%% 3.11
t_recov = linspace(0, 0.5, length(impulses_recov)+1); t_recov(end) = [];
f_recov = linspace(-fs_recov/2, fs_recov/2, length(t_recov)+1); f_recov(end) = [];
freqResponses = fft(impulses_recov, length(impulses_recov), 2);
for k = 1:3
    k = 2
    [text] = decoder(texts_recov(k, :));
    figure
    sgtitle(text)
    subplot(2, 1, 1)
    plot(t_recov, impulses_recov(k, :));
    title("Impulse Response")
    subplot(2, 1, 2)
    plot(f_recov, abs(fftshift(freqResponses(k, :)/fs_recov)));
end


%% 3.12
% do a voice recorder and copy path here, or use my beautiful voice 
%[voice, fs_voice] = audioread("C:\Users\turne\OneDrive\Documents\Sound recordings\Recording.m4a");
[voice, fs_voice] = audioread("voice.wav");
% convert to row vector and only use one channel
voice = voice(:, 1).';
voice = resample(voice, fs_recov, fs_voice);
% apply the impulse response as a convolution in the time domain
% 'same' will ensure the length is the same for both
voice_thing = conv(voice, impulses_recov(2, :), 'same');
sound(voice_thing, fs_recov)