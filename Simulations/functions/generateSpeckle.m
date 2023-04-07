function patternStack = generateSpeckle(nPixelX, nPixelY, pixelSize, NA, lambda, nImg)
% GENERATESPECKLE Summary of this function goes here
% Input parameters:
%   nPixelX: number of pixels in x; 
%   nPixelY: number of pixels in y;
%   pixelSize: pixel size in spatial domain
%   NA : NA of CTF;
%   lambda: wavelength of laser
%
% output parameters:
%   intensity_stack: speckle pattern

% nPixelX = 512;
% nPixelY = 512;
% pixelSize = 20; % unit:nm
% NA = 1.49;
% lambda = 488;
% nImg = 1000;
%   Detailed explanation goes here

F = @(x) fftshift(fft2(x));
iF = @(x) ifft2(ifftshift(x));

%% Generate the CTF
[Y,X]=meshgrid(1:nPixelY,1:nPixelX);
xc=floor(nPixelX/2+1);% the x-coordinate of the center
yc=floor(nPixelY/2+1);% the y-coordinate of the center
yr=Y-yc;
xr=X-xc;
R=sqrt((xr).^2+(yr).^2);% distance between the point (x,y) and center (xc,yc)

pixelNum=nPixelX;
rPixel=NA*pixelNum*pixelSize/lambda;
ctf=ones(pixelNum,pixelNum).*(R<=rPixel);
ctfSignificantPix=numel(find(abs(ctf)>eps(class(ctf))));
ifftscale=numel(ctf)/ctfSignificantPix;
apsf = iF(ctf);
ipsf = ifftscale*abs(apsf).^2;
OTF = real(F(ipsf));
OTF = OTF./max(abs(OTF(:)));

%% Pattern generation with random phase mask 
patternStack = zeros(nPixelX, nPixelY, nImg);

for iImg = 1:nImg

    % method 1
    random_map = rand(nPixelX, nPixelY);
    random_map(random_map<=0.95) = 0;
    random_map(random_map>0.95) = 1;
    random_mapf = F(random_map);
    pattern = abs(iF(random_mapf.*OTF));

    % method 2
%     random_mapf = exp(1j*rand(nPixelX, nPixelY)*100);
%     pattern = abs(iF(random_mapf.*ctf)).^2;
    
    patternStack(:,:,iImg) = pattern;
end
    
patternStack =  patternStack./max(patternStack(:));
figure;imagesc(pattern(:,:,1));colormap gray;axis square;title('pattern');colorbar;
figure;imagesc(mean(patternStack,3)); colormap gray;axis square;title('mean pattern');colorbar;
end


