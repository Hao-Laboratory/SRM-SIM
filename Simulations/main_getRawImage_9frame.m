% get raw image of star-like /particle sample��with noise��, illuminated by
% Sinusoidal fringes
% programmed by Shijie Tu
% 11830027@zju.edu.cn
% Dec 08, 2022

clear;close all;
addpath('.\functions\');
%% parameter input
lambdaEmi = 525; % Emission central wavelength, unit: nm.
NA = 1.49; % Numerical aperture of objective len.

load 'sample_star_like_80.mat';

[nPixelX, nPixelY] =  size(object);
pixelSize = 15.2;    % Pixel size in the imaging matrix (or 1-D array), unit: nm.

% parameters of illumination patterns
period = 76*4; % period of pattern, unit: nm.
nOri = 3; % number of pattern's orientations.
nPhase = 3; % number of pattern's phase.

% parameters of imaging environment
distortFactor = 0; % number for distortion
SNRdB_Sample = 15; % SNR, unit:dB,SNR = 10log10[var(image)/var(error image)]
SNRdB_Pattern = 50;
bkgLvl_Sample = 0.5;
bkgLvl_Pattern = 0.1;

fileNameSample = ['period-',num2str(period),' distortFactor-',num2str(distortFactor), ' SNR-', num2str(SNRdB_Sample),'dB bkg-',num2str(bkgLvl_Sample)];
fileNamePattern = ['period-',num2str(period),' distortFactor-',num2str(distortFactor), ' SNR-', num2str(SNRdB_Pattern),'dB bkg-',num2str(bkgLvl_Pattern)];
filePath = ['.\results\9 frame\',fileNameSample,'\']; % the saves path
mkdir(filePath);

save ([filePath,'simu_parameter.mat']);

%% initialize the sample matrix (2-D array)
sample = object./max(object(:));
imwrite(uint16(sample*65535), [filePath,'sample.tif']);

%% emission Point Spread Function(PSF)
[psf, OTF] = generatePSF(nPixelX,nPixelY,pixelSize,NA,lambdaEmi);
[psfBkg, ~] = generatePSF(nPixelX,nPixelY,pixelSize,0.5*NA,lambdaEmi);

%% results
% initialize the matrics to store the results
sampleImgs = zeros(nPixelX, nPixelY, nPhase*nOri);
patternImgs = zeros(nPixelX, nPixelY, nPhase*nOri);
excPatternImgs = zeros(nPixelX, nPixelY, nPhase*nOri);
% generate sinusoidal illumination pattern
excPattern = genernateSinuPattern(nPixelX, nPixelY, ...
    pixelSize, period, nOri, nPhase, distortFactor);

excWFTemp = distortImage(ones(4*nPixelX, 4*nPixelY), distortFactor);
excWF = excWFTemp(1.5*nPixelX+1:2.5*nPixelX, 1.5*nPixelY+1:2.5*nPixelY);

bkgSample = bkgLvl_Sample*conv2(excWF.*sample, psfBkg, 'same');
bkgPattern = bkgLvl_Pattern*conv2(excWF, psfBkg, 'same');

% calcualtion to get the raw images                                 
for iOri = 1 : nOri
    for iPhase = 1 : nPhase
        index = nPhase*(iOri-1)+iPhase;

        iExc = excPattern(:, :, iPhase, iOri);
        iSampleImg = conv2(iExc.*sample, psf, 'same');
        iSampleImg = addBkgNoise(iSampleImg, bkgSample, SNRdB_Sample);
        sampleImgs(:, :, index) = iSampleImg;  

        iPatternImg = conv2(iExc, psf, 'same');
        iPatternImg = addBkgNoise(iPatternImg, bkgPattern, SNRdB_Pattern);
        patternImgs(:, :, index) = iPatternImg;

        excPatternImgs(:, :, index) = iExc;
    end

end  % Get the image of the sample when the when the structed illumation pattern shifted.


imWriteStack(sampleImgs,[filePath,'sampleImgs ',fileNameSample,'.tif']);
imWriteStack(patternImgs,[filePath,'patternImgs ',fileNamePattern,'.tif']);
imWriteStack(excPatternImgs,[filePath,'excPatternImgs',' distortFactor-', num2str(distortFactor),'.tif']);

%% reconstruct image
SRMSIM = mean((sampleImgs-mean(sampleImgs,3)).*(excPatternImgs-mean(excPatternImgs,3)),3);

imWriteStack(SRMSIM, [filePath,'SRMSIM.tif']);
imWriteStack(mean(sampleImgs,3), [filePath,'WF.tif']);

