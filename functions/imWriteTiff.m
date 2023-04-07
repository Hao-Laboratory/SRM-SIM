function imWriteTiff(tifImage, fileName)
i = Tiff(fileName, 'w');
infostruct.ImageLength = size(tifImage, 1);
infostruct.ImageWidth = size(tifImage, 2);
infostruct.Photometric = Tiff.Photometric.MinIsBlack;
infostruct.BitsPerSample = 16;
infostruct.SampleFormat = Tiff.SampleFormat.UInt;
infostruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
i.setTag(infostruct)
i.write(uint16(65535*tifImage/max(tifImage(:))));
i.writeDirectory();
i.close();
end