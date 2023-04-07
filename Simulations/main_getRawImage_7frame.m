%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Spatially remodulated structured illumination microscopy    
%           Copyright (C) 2023 Shijie Tu  
%
% get raw image of star-like sample£¨with background and noise£©illuminated by
% sinusoidal patterns
% programmed by Shijie Tu
% 11830027@zju.edu.cn
%                                                   
% Please cite:                                      
%                                                   
% Shijie Tu et al., "High-speed spatially re-modulated structured
% illumination microscopy," DOI 10.1364/OL.485929                                    
%                                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;close all;
addpath('.\functions\');

% define fft, inverse fft and ihcoherent imaging function
F = @(x) fftshift(fft2(x));
iF = @(x) ifft2(ifftshift(x));
incoherentImg = @(object,OTF) real(iF(F(object).*OTF));

%% parameter input
lambdaEmi = 525; % Emission central wavelength, unit: nm.
NA = 1.49; % Numerical aperture of objective len.

load 'sample_star_like_80.mat';

[nPixelX, nPixelY] =  size(object);
pixelSize = 15.2;    % Pixel size in the imaging matrix (or 1-D array), unit: nm.

period = 76*4; % period of pattern, unit: nm.
nOri = 3; % number of pattern's orientations.
nPhase = 2; % number of pattern's phase.
distortFactor = 0; % number for distortion
SNRdB_Sample = 25; % SNR, unit:dB,SNR = 10log10[var(image)/var(error image)]
SNRdB_Pattern = 50;
bkgLvl_Sample = 0.5;
bkgLvl_Pattern = 0.1;

fileNameSample = ['period-',num2str(period),' distortFactor-',num2str(distortFactor), ' SNR-', num2str(SNRdB_Sample),'dB bkg-',num2str(bkgLvl_Sample)];
fileNamePattern = ['period-',num2str(period),' distortFactor-',num2str(distortFactor), ' SNR-', num2str(SNRdB_Pattern),'dB bkg-',num2str(bkgLvl_Pattern)];
filePath = ['.\results\7 frame\',fileNameSample,'\']; % the saves path
mkdir(filePath);

save ([filePath,'simu_parameter.mat']);

%% initialize the sample matrix (2-D array)
sample = object./max(object(:));
imwrite(uint16(sample*65535), [filePath,'sample.tif']);

%% emission Point Spread Function(PSF)
[psf, OTF] = generatePSF(nPixelX,nPixelY,pixelSize,NA,lambdaEmi);
[psfBkg, OTFBkg] = generatePSF(nPixelX,nPixelY,pixelSize,0.5*NA,lambdaEmi);

%% results
% initialize the matrics to store the results
sampleImgs = zeros(nPixelX, nPixelY, nPhase*nOri + 1);
patternImgs = zeros(nPixelX, nPixelY, nPhase*nOri + 1);
excPatternImgs = zeros(nPixelX, nPixelY, nPhase*nOri + 1);
% generate sinusoidal illumination pattern
excPattern = genernateSinuPattern(nPixelX, nPixelY, ...
    pixelSize, period, nOri, 2*nPhase, distortFactor, OTF);

excWFTemp = distortImage(ones(4*nPixelX, 4*nPixelY), distortFactor);
excWF = excWFTemp(1.5*nPixelX+1:2.5*nPixelX, 1.5*nPixelY+1:2.5*nPixelY);

bkgSample = bkgLvl_Sample.*incoherentImg(excWF.*sample, OTFBkg);
bkgPattern = bkgLvl_Pattern.*incoherentImg(excWF, OTFBkg);

WF = incoherentImg(excWF.*sample, OTF);
WF = addBkgNoise(WF, bkgSample, SNRdB_Sample);
WFPattern = incoherentImg(excWF, OTF);
WFPattern = addBkgNoise(WFPattern, bkgPattern, SNRdB_Pattern);


% calcualtion to get the raw images                                 
for iOri = 1 : nOri
    for iPhase = 1 : nPhase
        index = nPhase*(iOri-1)+iPhase;

        iExc = excPattern(:, :, iPhase, iOri);
        iSampleImg = incoherentImg(iExc.*sample, OTF);
        iSampleImg = addBkgNoise(iSampleImg, bkgSample, SNRdB_Sample);
        sampleImgs(:, :, index) = iSampleImg;  

        iPatternImg = incoherentImg(iExc, OTF);
        iPatternImg = addBkgNoise(iPatternImg, bkgPattern, SNRdB_Pattern);
        patternImgs(:, :, index) = iPatternImg;

        excPatternImgs(:, :, index) = iExc;
    end

end  % Get the image of the sample when the when the structed illumation pattern shifted.

sampleImgs(:,:,nPhase*nOri + 1) = WF;
patternImgs(:,:,nPhase*nOri + 1) = WFPattern;
excPatternImgs(:,:,nPhase*nOri + 1) = excWF;

imWriteStack(sampleImgs,[filePath,'sampleImgs ',fileNameSample,'.tif']);
imWriteStack(patternImgs,[filePath,'patternImgs ',fileNamePattern,'.tif']);
imWriteStack(excPatternImgs,[filePath,'excPatternImgs',' distortFactor-', num2str(distortFactor),'.tif']);

%% reconstruct image
SRMSIM = mean((sampleImgs-sampleImgs(:,:,end)).*(excPatternImgs-excPatternImgs(:,:,end)),3);

imWriteStack(SRMSIM, [filePath,'SRMSIM-initial.tif']);
imWriteStack(WF, [filePath,'WF.tif']);


