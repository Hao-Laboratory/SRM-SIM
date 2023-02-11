function distortedImage = distortImage(rawImage, k)

% details: https://ww2.mathworks.cn/help/symbolic/examples/developing-an-algorithm-for-undistorting-an-image.html
% https://blogs.mathworks.com/steve/2006/08/04/spatial-transformations-defining-and-applying-custom-transforms/

r = @(x) sqrt(x(:,1).^2 + x(:,2).^2);
% w = @(x) atan2(x(:,2), x(:,1));
% f = @(x) [sqrt(r(x)) .* cos(w(x)), sqrt(r(x)) .* sin(w(x))];
f = @(x) [x(:,1).*(1+k*r(x).^2), x(:,2).*(1+k*r(x).^2)];
g = @(x, unused) f(x);

tform = maketform('custom', 2, 2, [], g, []);
distortedImage = imtransform(rawImage, tform, 'UData', [-1 1], 'VData', [-1 1], ...
    'XData', [-1 1], 'YData', [-1 1]);