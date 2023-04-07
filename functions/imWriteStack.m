function imWriteStack(stack, filename)

im = Tiff(filename, 'w');

infostruct.ImageLength = size(stack, 1);
infostruct.ImageWidth = size(stack, 2);
infostruct.Photometric = Tiff.Photometric.MinIsBlack;
infostruct.BitsPerSample =16;
infostruct.SampleFormat = Tiff.SampleFormat.UInt;
infostruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
infostruct.Compression = Tiff.Compression.None;

maxValue = max(stack(:));
for k = 1:size(stack, 3)
    im.setTag(infostruct);
    img = stack(:,:,k);
    im.write(uint16(65535*img/maxValue));
    im.writeDirectory();
end

im.close();
end