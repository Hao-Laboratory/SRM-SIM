function excPattern = genernateSinuPattern(nPixelX, nPixelY, ...
    pixelSize, period, nOri, nPhase, distortFactor, OTF)

% GENERNATEPATTERN generate excitation sinusoidal pattern
% Input:
%   nPixelX: number of pixel in x direction
%   nPixelY: number of pixel in y direction
%   period:period of the sinusoidal pattern, unit:nm
%   nPatternOrientation: number of pattern's orientations
%   nPatternPhase: number of pattern's phase
%   OTF: optical transfer function
% Output:
%   excitationPattern: excitation pattern (2-D)
% Example:
%   genernatePattern(-500:500, -500:500, 200, 3, 3)

%%
F = @(x) fftshift(fft2(x));
iF = @(x) ifft2(ifftshift(x));

%%
xSize = 4*nPixelX; 
ySize = 4*nPixelY; 
xCoorArray = -floor(xSize / 2) : floor(xSize / 2 ) - 1 ; 
xCoorArray = xCoorArray * pixelSize; 
yCoorArray = -floor(ySize / 2) : floor(ySize / 2 ) - 1; 
yCoorArray = yCoorArray * pixelSize; 
[X,Y] = meshgrid(xCoorArray, yCoorArray);

PhaseDeg = period / nPhase.* (0 : nPhase -1 );
% patternPhaseDegrees = round(period / nPatternPhase).* [0 1.05 2]; % phase with error
OriDeg = 180/nOri*(0 : nOri -1 ); % Step size when the structed illumation pattern shifted. unit:nm.
% patternAngleDegrees = [0 66 120];% orientation with error
excPattern = zeros(nPixelX, nPixelY, nPhase, nOri);

for iOri = 1 : nOri    
    for iPhase = 1 : nPhase
        xShiftedCoordinateArray = X -  PhaseDeg(iPhase) * cosd(OriDeg(iOri));
        yShiftedCoordinateArray = Y -  PhaseDeg(iPhase) * sind(OriDeg(iOri));       
        xPeriod = period / cosd(OriDeg(iOri));
        yPeriod = period / sind(OriDeg(iOri));
        xPhaseArray = 2 * pi / xPeriod * xShiftedCoordinateArray;
        yPhaseArray = 2 * pi / yPeriod * yShiftedCoordinateArray;
        
        iExcitation = 1 + cos(xPhaseArray + yPhaseArray);
        
        % distort the pattern
        iExcitation = distortImage(iExcitation, distortFactor);

        excPattern(:, :, iPhase, iOri) = iExcitation(1.5*nPixelX+1:2.5*nPixelX, 1.5*nPixelY+1:2.5*nPixelY);
           
    end

end

excPatternf = F(excPattern);
excPattern = abs(iF(excPatternf.*OTF));

% figure;
% imagesc(xCoordinateArray, xCoordinateArray, mean(mean(excitationPattern,3),4));
% colorbar;
% title('mean of patterns')

end

