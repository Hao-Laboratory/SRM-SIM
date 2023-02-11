function outputImg = addBkgNoise(I, bkg, SNRdB)

% generate gaussian noise
varGauss = var(I(:))/10^(SNRdB/10);
gaussNoise = sqrt(varGauss)*randn(size(I));

% add bkg and noise to image
outputImg = I + bkg + gaussNoise;

end
