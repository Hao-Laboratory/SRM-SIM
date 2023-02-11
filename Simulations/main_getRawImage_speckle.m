% get raw image of star-like /particle sample£¨with noise£©, illuminated by
% speckle patterns
% programmed by Shijie Tu
% 11830027@zju.edu.cn
% Feb 4, 2023

clear; clc; close all
addpath('.\functions\');
%% parameter input
lambdaExc = 488; % Excitaion central wavelength, unit: nm.
lambdaEmi = 525; % Emission central wavelength, unit: nm.
NA = 1.49; % Numerical aperture of objective len.
pixelSize = 15.2; % Pixel size in the imaging matrix (or 1-D array), unit: nm.
nPattern = 30; % number of speckle pattern

load 'sample_star_like_80.mat';

[nPixelX, nPixelY] =  size(object);

% parameters of imaging environment
distortFactor = 0; % number for distortion
SNRdB_Sample = 15; % SNR, unit:dB,SNR = 10log10[var(image)/var(error image)]
SNRdB_Pattern = 50;
bkgLvl_Sample = 0.5;
bkgLvl_Pattern = 0.1;

filePath = ['.\results\speckle\', num2str(nPattern),' patterns\'];

fileNameSample = ['pattern-', num2str(nPattern),' distortFactor-',num2str(distortFactor), ' SNR-', num2str(SNRdB_Sample),'dB bkg-',num2str(bkgLvl_Sample)];
fileNamePattern = ['pattern-', num2str(nPattern),' distortFactor-',num2str(distortFactor), ' SNR-', num2str(SNRdB_Pattern),'dB bkg-',num2str(bkgLvl_Pattern)];
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
[psfBkg, ~] = generatePSF(nPixelX,nPixelY,pixelSize,0.5*NA,lambdaEmi);

%% results
% initialize the matrics to store the results
sampleImgs = zeros(nPixelX, nPixelY, nPattern);
patternImgs = zeros(nPixelX, nPixelY, nPattern);
excPatternImgs = zeros(nPixelX, nPixelY, nPattern);

% generate speckle illumination pattern
excPattern = generateSpeckle(nPixelX, nPixelY, pixelSize, NA, lambdaExc, nPattern);
excPattern = distortImage(excPattern , distortFactor);

excWF = ones(nPixelX, nPixelY);
excWF = distortImage(excWF, distortFactor);

bkgSample = bkgLvl_Sample*conv2(excWF.*sample, psfBkg, 'same');
bkgPattern = bkgLvl_Pattern*conv2(excWF, psfBkg, 'same');

% calcualtion to get the raw images                                 
for index = 1 : nPattern
        iExc = excPattern(:, :, index);
        iSampleImg = conv2(iExc.*sample, psf, 'same');
        iSampleImg = addBkgNoise(iSampleImg, bkgSample, SNRdB_Sample);
        sampleImgs(:, :, index) = iSampleImg;  

        iPatternImg = conv2(iExc, psf, 'same');
        iPatternImg = addBkgNoise(iPatternImg, bkgPattern, SNRdB_Pattern);
        patternImgs(:, :, index) = iPatternImg;

        excPatternImgs(:, :, index) = iExc;

end  % Get the image of the sample when the when the structed illumation pattern shifted.


imWriteStack(sampleImgs,[filePath,'sampleImgs ',fileNameSample,'.tif']);
imWriteStack(patternImgs,[filePath,'patternImgs ',fileNamePattern,'.tif']);
imWriteStack(excPatternImgs,[filePath,'excPatternImgs',' distortFactor-', num2str(distortFactor),'.tif']);

%% reconstruct image

SRMSIM = mean((sampleImgs-mean(sampleImgs,3)).*(excPatternImgs-mean(excPatternImgs,3)),3);

imWriteStack(SRMSIM, [filePath,'SRMSIM.tif']);
imWriteStack(mean(sampleImgs,3), [filePath,'WF.tif']);





    
