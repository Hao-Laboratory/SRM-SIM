function [ipsf, OTF] = generatePSFOTFGPU(xSize,ySize,pixelSize,NA,lambda)

[Y,X]=meshgrid(1:ySize,1:xSize);
xc=floor(xSize/2+1);% the x-coordinate of the center
yc=floor(ySize/2+1);% the y-coordinate of the center
yr=Y-yc;
xr=X-xc;
R=sqrt((xr).^2+(yr).^2);% distance between the point (x,y) and center (xc,yc)

%% Generate the PSF and OTF
pixelNum=xSize;
rPixel=NA*pixelNum*pixelSize/lambda;

ctf = ones(xSize,ySize,'gpuArray');
ctf = ctf.*(R<=rPixel);
ctfSignificantPix=numel(find(abs(ctf)>eps));
ifftscale=numel(ctf)/ctfSignificantPix;
apsf=fftshift(ifft2(ifftshift(ctf)));
ipsf=ifftscale*abs(apsf).^2;
OTF=real(fftshift(fft2(ifftshift(ipsf))));
OTF = OTF./max(abs(OTF(:)));

end
