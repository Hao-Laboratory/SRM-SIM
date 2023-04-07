%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Spatially remodulated structured illumination microscopy    
%           Copyright (C) 2023 Shijie Tu           
%                                                   
% Please cite:                                      
%                                                   
% Shijie Tu et al., "High-speed spatially re-modulated structured
% illumination microscopy," DOI 10.1364/OL.485929                                    
%                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ---------- clear workspace ---------- 
clc;clear;close all;

%% ---------- define FFT, inverse FFT, Winner Filter and fourier filtering ---------- 
F = @(x) fftshift(fft2(x));
iF = @(x) ifft2(ifftshift(x));
WnrFilter = @(OTF, factor) conj(OTF)./(abs(OTF).^2 + factor);
FourierFiltering = @(img, fourierFilter) real(iF(F(img).*fourierFilter));

%% ---------- data information ----------
addpath('.\functions\');
dataDir = '.\data\9-frame\';
sampleName = '488-mito-6pix.tif';
patternName = '488-pattern-6pix.tif';
saveLoc = [dataDir,'results-GPU\',sampleName(1:end-4),'\'];
mkdir(saveLoc);

M = 512; N = 512; % number of pixels in x of y
nPhase = 3; % number of phase shift in each direction: 3 for 9-frame SIM and 2 for 7-frame SIM
nDir = 3;
pixelSize = 65; % pixel size. unit: nm
lambdaEmi = 525; % fluorescence emission wavelength (emission maximum). unit: nm (670/525)
NA = 1.49;
  
%% ---------- reconstruction parameters setting ---------- 
scaleFactor = 2;
pixelSize = pixelSize/scaleFactor;
M = M*scaleFactor;
N = N*scaleFactor;

% [empirical value] 
ampPre = 1; % ranges from 0 to 1, depends on the background level of raw data, e.g. 1
sgimaPre = 1; % positive number，e.g. 1
ampPost = 0.75; % ranges from 0 to 1, depends on the background level of corvariance image, e.g. 0.75
sgimaPost = 0.01; % positive number e.g. 0.01
wnrFactorWF = 0.01; % wiener factor for pattern images pre-processing, positive number, e.g. 0.01
wnrFactorSIM = 0.001; % wiener factor for SRM-SIM post-processing, positive number, e.g. 0.001
shadingFactor = 0.001;% wiener factor for shading correction, positive number, e.g. 0.001

% check parameters input
if (sgimaPre<=0 || sgimaPost<=0 || wnrFactorWF<=0 || wnrFactorSIM<=0 ...
        || shadingFactor<=0 || ampPre<0 || ampPost<0)
    error("Wrong parameters input: found negative value!!!")
end

save([saveLoc,'SRMSIM parameters.mat']);

%% ---------- preparations (can be reused for consecutive reconstructions) ---------- 

% generate OTF and wiener filger for WF
[~, OTFWF] = generatePSFOTFGPU(M, N, pixelSize, NA, lambdaEmi);
wnrFilterWF = WnrFilter(OTFWF, wnrFactorWF);

% load pattern for SRM-SIM
pattern = imresize(single(imReadStackGPU([dataDir,patternName])),scaleFactor,'Method','bilinear');
nFrames = size(pattern,3);

% check data input
if (nPhase==3 && nFrames~=9) || (nPhase==2 && nFrames~=7)
    error("Wrong nPhase input!!!")
end

% pre-processing on pattern images
for i = 1 : nFrames
    pattern(:,:,i) = FourierFiltering(pattern(:,:,i), wnrFilterWF);
end

% calculate (pattern - WFPattern)
if nPhase == 3
    patternSubMean = pattern - mean(pattern,3);
else
    patternSubMean = pattern - pattern(:,:,end);
end

% get shading filter (optional and useful for distortion correction)
shading = mean(patternSubMean.*patternSubMean,3);
shading = shading./max(shading(:));
shadingFilter = WnrFilter(shading, shadingFactor); 

% get attenuation function used in pre/post-processing
[Y,X] = ndgrid(1:M,1:N);
rad = gpuArray(single(hypot(X-floor(N/2+1),Y-floor(M/2+1))));
cyclesPerMicron = 1/(N*pixelSize/1000);
cycl = rad.*cyclesPerMicron;
attFunPre = 1-ampPre*exp(-cycl.^2./(2*sgimaPre^2));
attFunPost = 1 - ampPost*exp(-cycl.^2./(2*(sgimaPost)^2));

% get notch filter used in pre-processing
notchFilter = attFunPre.*conj(OTFWF);
notchFilter = notchFilter./max(notchFilter(:));

% get effective OTF for post-processing (final Wiener deconvolution)
notchedEffOTF = OTFWF.*notchFilter; % effective OTF after pre-processing (notch filtering)
notchedEffOTF = notchedEffOTF./max(notchedEffOTF(:));
shiftedEffOTF = gpuArray(single(zeros(size(notchedEffOTF))));
kr = gpuArray(single(zeros(nDir,1)));
for iDir = 1:nDir
    [kx, ky] = getVectorPerDir(pattern(:,:,nPhase*(iDir-1)+1), OTFWF); % get frequency vector per orientation
    kr(iDir) = sqrt(kx.^2 + ky.^2);
    OTFtemp = circshift(notchedEffOTF, [-kx -ky]) + circshift(notchedEffOTF, [kx ky]);
    shiftedEffOTF = shiftedEffOTF + OTFtemp;
end
shiftedEffOTF = shiftedEffOTF./max(shiftedEffOTF(:));

% get apodization function for SRM-SIM
kcWF = 2*NA*M*pixelSize/lambdaEmi;
krPattern = mean(kr(:));
kcSIM = kcWF + krPattern;
apoFunSIM = gpuArray(single(cos(pi*rad/(2*kcSIM)))).*attFunPost;  
apoFunSIM( rad > kcSIM ) = 0; 
apoFunSIM = apoFunSIM./max(apoFunSIM(:));

% get filter used in post-processing
wnrFilterSIM = apoFunSIM.*WnrFilter(shiftedEffOTF, wnrFactorSIM);

%% ---------- load raw SIM images and create the  pseudo-widefield image ---------- 
sample = imresize(single(imReadStackGPU([dataDir,sampleName])),scaleFactor,'Method','bilinear');

if nPhase == 3
    WF = mean(sample,3);
else 
    WF = sample(:,:,end);
end

%% ---------- main reconstruction procedure of SRM-SIM ---------- 
tic;
% pre-processing on sample images
for i = 1 : nFrames
    sample(:,:,i) = real(ifft2(ifftshift((fftshift(fft2(sample(:,:,i)))).*notchFilter)));
end

% get super-resolutoin image with spatial remodulation
if nPhase == 3
    sampleSubMean = sample - mean(sample,3); % calculate (sample - WF)
else
    sampleSubMean = sample - sample(:,:,end);
end
SRMImg = mean(sampleSubMean.*patternSubMean,3).*shadingFilter;

% post-processing (final Wiener deconvolution)
SRMSIM = real(ifft2(ifftshift((fftshift(fft2(SRMImg))).*wnrFilterSIM)));
toc;

%% ---------- save the final super-resolved image ---------- 
imWriteTiff(gather(WF),[saveLoc,'WF.tif']);
imWriteTiff(gather(SRMSIM),[saveLoc,'SRMSIM.tif']);
