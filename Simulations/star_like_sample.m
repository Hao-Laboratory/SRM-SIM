%%
nPixel = 512;
nBlockPixel = 20; % the vale of the block pixel is zero
x = -nPixel/2 : nPixel/2-1;
y = x;
[X, Y] = meshgrid(x, y);
theta = atan2(Y, X);
object = 1 + cos(80*theta);
object(:, 1:nBlockPixel) = 0;
object(:, nPixel-nBlockPixel+1 : nPixel) = 0;
object(1:nBlockPixel, :) = 0;
object(nPixel-nBlockPixel+1 : nPixel, :) = 0;

figure;
imagesc(object);axis square;
colorbar;

save sample_star_like_80 object;