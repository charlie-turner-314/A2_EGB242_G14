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

%% 3.1 - Spectrum Analysis
Ts = 1/fs;                                                  % Sampling Period [s]
samples = length(muxSignal);                                % Number of samples
t = linspace(0, samples/fs, samples+1); t(end) = [];        % Time vector
MUXSIG = fft(muxSignal);                                    % Fourier Transform
f = linspace(-fs/2, fs/2, samples+1); f(end)=[];            % Frequency Vector

% Plot signal against time
figure('Name',"Recieved Signal"), hold on;
subplot(1, 2, 1)
plot(t, muxSignal);                                         % Time Domain
xlabel("Time [s]"), ylabel("Amplitude"), title("Time Domain of muxSignal")
subplot(1, 2, 2)
plot(f, abs(fftshift(MUXSIG/fs)));                          % Frequency Domain
xlabel("Frequency [Hz]"), ylabel("Magnitude"), title("Magnitude Spectrum of muxSignal")

%% 3.2 - Multiplexing Parameters
% We can note that there are 6 clear frequency shifts from above
MUXSIGshift = fftshift(MUXSIG);                             % Shifted Fourier Transform to correspond to f
% Find the 6 frequency shifts, considering only positive freqs
[~, fshift] = findpeaks(abs(MUXSIGshift(f>0)), f(f>0), 'NPeaks',6, 'SortStr','descend');
fshift = sort(fshift);                                      % Sort in asc order
fshift_indices = find(ismember(f, fshift));                 % Find incices in f
Mag = abs(MUXSIGshift(fshift_indices));                     % Corresponding Magnitudes
Phase = angle(MUXSIGshift(fshift_indices));                 % Corresponding Phase [rads]    
clear MUXSIGshift;     % Done with this
%% 3.3 - Unshift Each Signal &  3.4 - Compute Fourier Transform
xdm = FDMDemux(muxSignal, t, Mag, fshift, Phase);           % Demultiplex using known values
XDM = fft(xdm, length(xdm), 2);                             % Row-wise FT

figure('Name', "Demultiplexed Streams"), hold on;
streams = size(xdm, 1);
cols = 3;
for k = 1:streams                                           % Plot each stream, 3 colums, 2 rows per stream
    subplot(streams/3 * 2, 3, k + mod(k-1, streams) - mod(k-1, 3)) % Don't mind the stupid way of laying out subpots
    plot(t, xdm(k, :));                                     % Time domain
    xlabel("Time [s]"), ylabel("Amplitude"), title(sprintf("Time Domain - Stream %i", k))
    subplot(streams/3 * 2, 3, k + mod(k-1, streams) - mod(k-1, 3) + 3)
    plot(f, abs(fftshift(XDM(k, :)/fs)));                   % Frequency Domain
    xlabel("Frequency [Hz]"), ylabel("Magnitude"), title(sprintf("Magnitude Spectrum of Stream %i", k))
end
%% 3.4, 3.5 - Analyse TFs
factorTF(sys(1))
% shows its  a time multiplexed exponential -> (1/(s+a)^2)
% No zeros, poles are real and neg so looks good for exponential decay
factorTF(sys(2))
% Neg poles, larger than sys 1 so may have a quicker decay
% has a zero at the origin, indicates low frequency attenuation
factorTF(sys(3))
% Complex poles, can infer oscillation within impulse response
% zeros at the origin, high pass filter 
factorTF(sys(4))
% Positive poles, unstable characteristics
%% 3.6 - Inspect Filters
% Calculate frequency responses of systems at our frequencies of interest
freqResponse = zeros([length(sys), length(f)]);             % Preallocate for performance
for k = 1:length(sys)                                       % Evaluate the freq response of each system
    freqResponse(k, :)= freqresp(sys(k), f, 'Hz');          % Substitute frequencies for S and evaluate output
end
figure('Name', "System Analysis", "Position", [50 50 1400 600]), hold on
subplot(2, 2, 1)
pzmap(sys(1), sys(2), sys(3), sys(4))
legend("Sys 1", "Sys 2", "Sys 3", "Sys 4")
% Clearly purple (sys4) is undesirable due to positive poles
subplot(2, 2, 2)
impulse(sys(1), sys(2), sys(3))
% Yellow and red go neg which doesn't look ideal
subplot(2, 2, 3), hold on
for k = 1:size(freqResponse, 1)
    plot(f, real(freqResponse(k, :)))
end
% Sys 4 clearly not ideal, sys 1 does what we want
xlabel("Frequency"), ylabel("Magnitude"), title("Frequency Response in f")
subplot(2, 2, 4)
bodemag(sys(1), sys(2), sys(3), sys(4))
% confirms 1 is a low pass
%% 3.7 - Choosing A Filter
% Clearly (sys(1)) is best, as it is a low pass with a close to ideal
% step and impulse response
%% 3.8  Apply the filter
syms s;
[Num,Den] = tfdata(sys(1),'v');                             % Extract numerator and denominator
sys1_TF = poly2sym(Num, s)/poly2sym(Den, s);                % Get symbolic transfer function
sys1_IR = matlabFunction(ilaplace(sys1_TF));                % Usable function for the impulse response
impresp = sys1_IR(t);                                       % Discrete impulse response vector
freqRes = fft(impresp);                                     % Discrete Frequency response
MSG = zeros(size(XDM));                                     % Inilialise size For Performance
for k = 1:size(XDM, 1)
    MSG(k, :) = freqRes .* XDM(k, :);                       % Apply the Transfer Function
end
msg = ifft(MSG, length(MSG), 2);                            % Convert to time domain
for k = 1:size(msg,1)
    msg(k, :) = msg(k, :) - mean(msg(k,:));                 % remove dc bias (assuming mean of signals is 0);
    msg(k, :) = msg(k, :) / (max(abs(msg(k, :))));          % Scale to [-1 1]
end

% Plot the new streams
figure('Name', "Time and Freq Domain Of Filtered and Rescaled Streams"), hold on     
for k = 1:streams
    subplot(streams/3 * 2, 3, k + mod(k-1, streams) - mod(k-1, 3))
    plot(t, msg(k, :))
    title(sprintf("Time Domain of Stream %i", k)), xlabel("Time [s]"), ylabel("Amplitude")
    ylim([-1, 1])
    subplot(streams/3 * 2, 3, k + mod(k-1, streams) - mod(k-1, 3) + 3)
    plot(f, real(abs(fftshift(MSG(k, :)/fs))))
    title(sprintf("Magnitude Spectrum of Stream %i", k)), xlabel("Frequency [Hz]"), ylabel("Magnitude")
end
%% 3.9
FS_recov = 16e3;                                                % ADC sample rate [Fs]
impulses_recov = resample(msg((1:2:end), :).', FS_recov, fs).'; % Resample each impulse response
%% 3.10
texts_recov = msg((2:2:end), :);
%% 3.11
t_recov = linspace(0, 0.5, length(impulses_recov)+1); t_recov(end) = [];
f_recov = linspace(-FS_recov/2, FS_recov/2, length(t_recov)+1); f_recov(end) = [];
IMPULSES_RECOV = fft(impulses_recov, length(impulses_recov), 2);
for k = 1:3
    [text] = decoder(texts_recov(k, :));
    figure
    sgtitle(text)
    subplot(1, 2, 1)
    plot(t_recov, impulses_recov(k, :));
    ylim([-1 1]), xlim([0 0.3])
    title("Impulse Response"), xlabel("Time [s]"), ylabel("Amplitude")
    subplot(1, 2, 2)
    plot(f_recov, abs(fftshift(IMPULSES_RECOV(k, :)/FS_recov)));
    title("Frequency Response"), xlabel("Frequency [Hz]"), ylabel("Magnitude")
end


%% 3.12
% do a voice recorder and copy path here, or use my beautiful voice 
%[voice, fs_voice] = audioread("C:\Users\turne\OneDrive\Documents\Sound recordings\Recording.m4a");
[voice, fs_voice] = audioread("voice.wav");
% convert to row vector and only use one channel
voice = voice(:, 1).';
voice = voice / (max(abs(voice)));
voice = resample(voice, FS_recov, fs_voice);
% apply the impulse response as a convolution in the time domain
voiceTest = zeros([size(impulses_recov, 1) length(voice)]);
for k=1:size(impulses_recov, 1)
    c = conv(voice, impulses_recov(k, :), 'same');
    voiceTest(k, :) = c / max(abs(c(:)));
end

t_voice = linspace(0, length(voice)/fs, length(voice)+1); t_voice(end) = [];
f_voice = FS_recov/2 * linspace(-1, 1, length(t_voice)+1); f_voice(end)=[];
figure('Name', 'Voice Testing')
subplot(4, 2, 1)
plot(t_voice, voice);
title("Time Domain (Voice recording)", 'Interpreter','latex'), xlabel("Time [s]"), ylabel("Amplitude")
subplot(4, 2, 2)
plot(f_voice, abs(fftshift(fft(voice)/FS_recov)))
title("Magnitude Spectrum (Voice Recording)", 'Interpreter','latex'), xlabel("Frequency [Hz]"), ylabel("Magnitude")
rooms = ["laboratory", "livingroom", "kitchen"];
for k=1:3
    row = k;
    subplot(4, 2, row*2+1)
    plot(t_voice, voiceTest(k, :))
    title(sprintf("Time Domain (Voice recording $*$ %s impulse response)", rooms(k)), 'Interpreter','latex'), xlabel("Time [s]"), ylabel("Amplitude")
    subplot(4, 2, row*2+2)
    plot(f_voice, abs(fftshift(fft(voiceTest(k, :))/FS_recov)))
    title(sprintf("Magnitude Spectrum (Voice recording $\\times$ %s frequency response)", rooms(k)), 'Interpreter','latex'), xlabel("Frequency [Hz]"), ylabel("Magnitude")
end

% sound(voiceTest(3, :), FS_recov)
% sound(voice, FS_recov)
