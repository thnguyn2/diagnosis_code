function res=extract_histogram_from_a_single_core_g3_vs_nm(curfilename,textonhistfilename,roicoordfilename,param)
%Inputs:
%   curfilename: name of the current label image
%   textonhistfilename: name of the histogram of texton indices
%   param: parameters for refining the gland and luminal regions
%Outputs:
%        res.nlumen_arr: an array of number of lumen area for each gland in
%        the region
%        res.distortion_arr: an array of distortion ratio for each gland.
%        It is the ratio of the gland's perimeter over pi*gland equi radius
%        res.area_arr:  an array of gland areas
%        res.nglands:   number of glands
%        res.mean_nlum: mean number of lumen;
%        res.mean_dist: mean distortion of the gland

        tempdata= load(textonhistfilename,'histim');
        histim = tempdata.histim;
        tempdata = load(roicoordfilename);
        r1 = tempdata.r1;r2 = tempdata.r2;c1=tempdata.c1;c2 = tempdata.c2;
        histim = histim(r1:r2,c1:c2,:);
        nrows = r2-r1+1;
        ncols = c2-c1+1;
        [x,y]=meshgrid([c1:c2],[r1:r2]);   
        lblim = imread(curfilename);
        glandim = (lblim==1);          
        glandim = imfill(glandim,'hole');
        stromaim = (lblim==2);
        %Extracting areas around the glands and look at the histogram
        glandstrand2=imdilate(glandim,strel('disk',param.basalwidth))-imerode(glandim,strel('disk',param.basalwidth));
        %glandstrand2 = imdilate(stromaim,strel('disk',param.basalwidth)).*glandstrand2;
        basalidx = find(glandstrand2==1);
        res.c1 = c1;
        res.c2 = c2;
        res.r1 = r1;
        res.r2 = r2;
        res.nrows = nrows;
        res.ncols = ncols;
        res.basalidx = basalidx; %Coordinate of the data
        res.glandidx = find(glandim==1);
        ntextons = size(histim,3);
        basalhistdata = zeros(length(basalidx),ntextons); %Get the 3D data where each row is one histogram at 1 pixel
        glandhistdata = zeros(length(res.glandidx),ntextons); %Get the 3D data where each row is one histogram at 1 pixel
        for bandidx = 1:ntextons
            curband = histim(:,:,bandidx);
            basalhistdata(:,bandidx)=curband(basalidx);
            glandhistdata(:,bandidx) = curband(res.glandidx);
        end
        res.basalhistdata = basalhistdata;
        res.glandhistdata = glandhistdata;
         
end
