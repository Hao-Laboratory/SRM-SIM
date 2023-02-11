function ImStack = imReadStack(img)
% This function can import 3D image stack.
% The bit depth is limited to 8,16,32-bit

imgInfo = imfinfo(img);
imgRow = imgInfo(1).Height;
imgCol = imgInfo(1).Width;
imgDepth = length(imgInfo);
imgBitDepth = ['uint',num2str(imgInfo(1).BitDepth)];
ImStack = zeros(imgRow, imgCol, imgDepth,imgBitDepth);
for ii = 1 : imgDepth
    ImStack(:,:,ii) = imread(img, ii);
end

end