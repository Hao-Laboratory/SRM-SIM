function patternStack = generateMultifocal(nPixelX, nPixelY, pixelSize, NA, lambda, dotPeriodPx, windowSizePx, N_shiftx, N_shifty, shiftStepPx)
% GENERATEMULTIFOCAL Summary of this function goes here
% Input parameters:
%   nPixelX: number of pixels in x; 
%   nPixelY: number of pixels in y;
%   pixelSize: pixel size in spatial domain
%   NA : NA of CTF;
%   lambda: wavelength of laser
%
% output parameters:
%   intensity_stack: Multifocal pattern
% clear;
% nPixelX = 512;
% nPixelY = 512;
% pixelSize = 20; % unit:nm
% NA = 1.49;
% lambda = 488;
% dotPeriodPx = 24;
% windowSizePx = 12;
% N_shiftx = 3;
% N_shifty = 3;
% shiftStepPx = 8;
%   Detailed explanation goes here

F = @(x) fftshift(fft2(ifftshift(x)));
iF = @(x) fftshift(ifft2(ifftshift(x)));

%% Generate the OTF
xSize = 4*nPixelX;
ySize = 4*nPixelY;
[Y,X]=meshgrid(1:ySize,1:xSize);
xc=floor(xSize/2+1);% the x-coordinate of the center
yc=floor(ySize/2+1);% the y-coordinate of the center
yr=Y-yc;
xr=X-xc;
R=sqrt((xr).^2+(yr).^2);% distance between the point (x,y) and center (xc,yc)

pixelNum=xSize;
rPixel=NA*pixelNum*pixelSize/lambda;
ctf=ones(pixelNum,pixelNum).*(R<=rPixel);
ctfSignificantPix=numel(find(abs(ctf)>eps(class(ctf))));
ifftscale=numel(ctf)/ctfSignificantPix;
apsf = iF(ctf);
ipsf = ifftscale*abs(apsf).^2;
OTF = real(F(ipsf));
OTF = OTF./max(abs(OTF(:)));

%% Pattern generation with random phase mask 
%
pattern = 0.01*ones(ySize,xSize);

num_period_y = floor(ySize/dotPeriodPx);
num_period_x = floor(xSize/dotPeriodPx);

for i = 1:num_period_y
    for j = 1:num_period_x
        idx_y = i*dotPeriodPx - dotPeriodPx/2;
        idx_x = j*dotPeriodPx - dotPeriodPx/2;
        pattern((idx_y-(windowSizePx/2-1)):(idx_y+(windowSizePx/2)),(idx_x-(windowSizePx/2-1)):(idx_x+(windowSizePx/2))) = 1;
    end
end

figure;imagesc(pattern);colormap gray;axis image;axis off;

%
patternf = F(pattern);
pattern = abs(iF(patternf.*OTF));

% figure;imagesc(log(1+abs(patternf)));colormap gray;axis image;axis off;
% figure;imagesc(pattern);colormap gray;axis image;axis off;


%% Generate simulated data with angular/pixel shift
% Pixel shift part

nImg = N_shiftx*N_shifty;

pixel_shiftx = (-(N_shiftx-1)/2:(N_shiftx-1)/2).*shiftStepPx;
pixel_shifty = (-(N_shifty-1)/2:(N_shifty-1)/2).*shiftStepPx;

[pixel_shiftyy,pixel_shiftxx] = meshgrid(pixel_shifty,pixel_shiftx);

pixel_shift_stack = zeros(2,nImg);

pixel_shift_stack(1,:) = pixel_shiftyy(:);
pixel_shift_stack(2,:) = pixel_shiftxx(:);

patternStack = zeros(nPixelY, nPixelX, nImg);

for i = 1:nImg
    
    % Pixel shift part
    patternStack(:,:,i) = pattern(1.5*nPixelY+1+pixel_shift_stack(1,i):2.5*nPixelY+pixel_shift_stack(1,i), ...
        1.5*nPixelX+1+pixel_shift_stack(2,i):2.5*nPixelX+pixel_shift_stack(2,i));

end

patternStack =  patternStack./max(patternStack(:));
figure;imagesc(mean(patternStack,3)); colormap gray;axis square;title('mean');colorbar;
end


