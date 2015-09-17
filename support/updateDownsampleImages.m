function updateDownsampleImages(tifflist)
    %Update downsample images
    %Input: tifflist: a list of downsample images of the tiff files
    nfiles = size(tifflist,1);
    for fileidx=1:nfiles
        disp(['Downsampling: ' tifflist{fileidx}]);
        fullsizeimname = tifflist{fileidx};
        resizeimname=strcat(fullsizeimname(1:end-4),'_small.tif');
        if (~exist(resizeimname,'file'))
                im = imread(fullsizeimname);
                im = cast(im,'uint16');
                im = imresize(im,[3072 3072]);
                figure(1);
                imagesc(im);
                writeTIFF(im,resizeimname,'uint16');
                disp('Saving...');
                disp('Done.');
        end    
    end
end