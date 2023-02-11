function [kx, ky] = getVectorPerDir(rawData, OTF)
% Estimates the SIM illumination pattern in one orientation by peak search

% normalise raw data
rawData = rawData./max(rawData(:));

rawDataFFT = fftshift(fft2(rawData));

% mask 
mask = 1-OTF;
firstOrder = log(1+abs(rawDataFFT.*mask));

[~, index] = sort(firstOrder(:));
[yPosFirst, xPosFirst] = ind2sub(size(firstOrder),index(end-1:end));
kx = round((xPosFirst(2) - xPosFirst(1))/2);
ky = round((yPosFirst(2) - yPosFirst(1))/2);

disp(['frequency-vector = (' num2str(kx) ',' num2str(ky) ')']);

end


