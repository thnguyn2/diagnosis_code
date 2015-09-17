function writeTIFF(data, filename,type)
% writeTIFF(data, filename)
% writes data as a multi-channel TIFF with single prec. float pixels
   t = Tiff(filename, 'w');
   tagstruct.ImageLength = size(data, 1);
   tagstruct.ImageWidth = size(data, 2);
   tagstruct.Compression = Tiff.Compression.None;
   %tagstruct.Compression = Tiff.Compression.LZW;        % compressed
   tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
   
   if (nargin==2)
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        tagstruct.BitsPerSample =  32;                        % float data
   elseif (strcmp(type,'uint16'))
         tagstruct.SampleFormat = Tiff.SampleFormat.UInt; %See tiff class
         tagstruct.BitsPerSample =  16;
   elseif (strcmp(type,'uint8'))       
        tagstruct.SampleFormat = Tiff.SampleFormat.Int; %See tiff class
         tagstruct.BitsPerSample =  8;
   end
   tagstruct.SamplesPerPixel = size(data,3);
   tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
   t.setTag(tagstruct);
   if (nargin==2)
       t.write(single(data));
   elseif (strcmp(type,'uint16'))
       t.write(uint16(data));
   elseif (strcmp(type,'int8'))
       t.write(int8(data));
   t.close();
end