%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Spatially remodulated structured illumination microscopy    
%           Copyright (C) 2023 Shijie Tu  
%
% get raw image of star-like sample£¨with background and noise£©illuminated by
% multifocal patterns
% programmed by Shijie Tu
% 11830027@zju.edu.cn
%                                                   
% Please cite:                                      
%                                                   
% Shijie Tu et al., "High-speed spatially re-modulated structured
% illumination microscopy," DOI 10.1364/OL.485929                                    
%                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all
addpath('.\functions\');

% define fft, inverse fft and ihcoherent imaging function
F = @(x) fftshift(fft2(x));
iF = @(x) ifft2(ifftshift(x));
incoherentImg = @(object,OTF) real(iF(F(object).*OTF));

%% parameter input
lambdaExc = 488; % Excitaion central wavelength, unit: nm.
lambdaEmi = 525; % Emission central wavelength, unit: nm.
NA = 1.49; % Numerical aperture of objective len.
pixelSize = 15.2;    % Pixel size in the imaging matrix (or 1-D array), unit: nm.
dotPeriodPx = round(76*3/pixelSize); % period of dot, unit: pixel
windowSizePx = round(76/pixelSize);% size of dot, unit: pixel
N_shiftx = 3; % number of shift step in x direction
N_shifty = 3; % number of shift step in y direction
shiftStepPx = dotPeriodPx./N_shiftx; % shift step size, unit: pixel

load 'sample_star_like_80.mat';

[nPixelX, nPixelY] =  size(object);

% parameters to get speckle pattern
nPattern = N_shiftx.*N_shifty; % number of multifocal pattern

% parameters of imaging environment
distortFactor = 0; % number for distortion
SNRdB_Sample = 25; % SNR, unit:dB,SNR = 10log10[var(image)/var(error image)]
SNRdB_Pattern = 50;
bkgLvl_Sample = 0.5;
bkgLvl_Pattern = 0.1;

filePath = ['.\results\multifocal\', num2str(nPattern),' patterns\'];

fileNameSample = ['multifocal-pattern-', num2str(nPattern),' distortFactor-',num2str(distortFactor), ' SNR-', num2str(SNRdB_Sample),'dB bkg-',num2str(bkgLvl_Sample)];
fileNamePattern = ['multifocal-pattern-', num2str(nPattern),' distortFactor-',num2str(distortFactor), ' SNR-', num2str(SNRdB_Pattern),'dB bkg-',num2str(bkgLvl_Pattern)];
mkdir(filePath);

save ([filePath,'simu_parameter.mat']);

%% values of key parameters based on the inputs
xCoorArray = -floor(nPixelX / 2) : floor(nPixelX / 2 ) - 1 ; 
xCoorArray = xCoorArray * pixelSize; 
yCoorArray = -floor(nPixelY / 2) : floor(nPixelY / 2 ) - 1; 
yCoorArray = yCoorArray * pixelSize; 

[X,Y] = meshgrid(xCoorArray, yCoorArray);

%% initialize the sample matrix (2-D array)
sample = object./max(object(:));
imwrite(uint16(sample*65535), [filePath,'sample.tif']);

%% emission Point Spread Function(PSF)
[psf, OTF] = generatePSF(nPixelX,nPixelY,pixelSize,NA,lambdaEmi);
[psfBkg, OTFBkg] = generatePSF(nPixelX,nPixelY,pixelSize,0.5*NA,lambdaEmi);

%% results
% initialize the matrics to store the results
sampleImgs = zeros(nPixelX, nPixelY, nPattern);
patternImgs = zeros(nPixelX, nPixelY, nPattern);
excPatternImgs = zeros(nPixelX, nPixelY, nPattern);

% generate speckle illumination pattern
excPattern = generateMultifocal(nPixelX, nPixelY, pixelSize, NA, lambdaExc, ...
    dotPeriodPx, windowSizePx, N_shiftx, N_shifty, shiftStepPx);
excPattern = distortImage(excPattern , distortFactor);

excWF = ones(nPixelX, nPixelY);
excWF = distortImage(excWF, distortFactor);

bkgSample = bkgLvl_Sample.*incoherentImg(excWF.*sample, OTFBkg);
bkgPattern = bkgLvl_Pattern.*incoherentImg(excWF, OTFBkg);

% calcualtion to get the raw images                                 
for index = 1 : nPattern
        iExc = excPattern(:, :, index);
        iSampleImg = incoherentImg(iExc.*sample, OTF);
        iSampleImg = addBkgNoise(iSampleImg, bkgSample, SNRdB_Sample);
        sampleImgs(:, :, index) = iSampleImg;  

        iPatternImg = incoherentImg(iExc, OTF);
        iPatternImg = addBkgNoise(iPatternImg, bkgPattern, SNRdB_Pattern);
        patternImgs(:, :, index) = iPatternImg;

        excPatternImgs(:, :, index) = iExc;

end  % Get the image of the sample when the when the structed illumation pattern shifted.


imWriteStack(sampleImgs,[filePath,'sampleImgs ',fileNameSample,'.tif']);
imWriteStack(patternImgs,[filePath,'patternImgs ',fileNamePattern,'.tif']);
imWriteStack(excPatternImgs,[filePath,'excPatternImgs',' distortFactor-', num2str(distortFactor),'.tif']);

%% reconstruct image
WF = mean(sampleImgs,3);
SRMSIM = mean((sampleImgs-mean(sampleImgs,3)).*(excPatternImgs-mean(excPatternImgs,3)),3);

imWriteStack(SRMSIM, [filePath,'SRMSIM-initial.tif']);
imWriteStack(WF, [filePath,'WF.tif']);





    
