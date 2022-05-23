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
