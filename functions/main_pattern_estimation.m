
clc;clear;close all;

% SIM data
SIM_path = '..\data\9-frame\'; % path to data files
SIM_name = '488-mito-6pix'; % sim data 
SIM_suffix = '.tif';
HEIGHT = 512;
WIDTH = 512;

pixelSize = 65; % pixel size. unit: nm
lambdaEmi = 525; % fluorescence emission wavelength (emission maximum). unit: nm (670/525)
NA = 1.49;

data = double(imReadStack([SIM_path SIM_name SIM_suffix]));

% generate OTF
[~, otf] = generatePSFOTF(WIDTH, HEIGHT, pixelSize, NA, lambdaEmi);
otf = otf/max(otf(:));

% Pattern estimation
reg_beta = 1e-3;
reg_delta = 1e-3;
max_itr = 50;

tic;
[Ip_est,I_mean,I_mean_dec] = IPE(data, otf, max_itr, reg_beta, reg_delta);
toc;

% imWriteStack(I_mean,[SIM_path SIM_name '-WF' SIM_suffix]);
% imWriteStack(I_mean_dec,[SIM_path SIM_name '-WF-dec' SIM_suffix]);
imWriteStack(Ip_est,[SIM_path SIM_name '-IPE' SIM_suffix]);