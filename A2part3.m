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
%% Setup
B = 8e3;                    % Bandwidth  [Hz]

%% 3.1
Ts = 1/fs;                  % Sampling Period
t = 0:Ts:(length(muxSignal)/fs); t(end) = [];
% Plot t
figure(1), hold on;
subplot(1, 2, 1)
plot(t, muxSignal);
xlabel("Time [s]"), ylabel("Amplitude"), title("Time Domain of muxSignal")

MUXSIG = fft(muxSignal);
f = linspace(-fs/2, fs/2, length(t)+1); f(end)=[];
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
figure(2), hold on;
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
end
%% 3.5
% bode plot shows low pass is number 1. Make some stuff ap about the
% denominator lol
% factorTF(sys(1)) % shows its  a time multiplexed exponential -> (1/(s+a)^2)
% figure(), title("System 1 Freq Response")
% h = bodeplot(sys(4));
% setoptions(h, 'FreqUnits', 'Hz')
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
%% 3.9
%% 3.10
%% 3.11
%% 3.12