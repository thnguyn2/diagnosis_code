function res=process_single_core_g3_vs_nm(curfilename,param)
%Inputs:
%   curfilename: name of the current label image
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

        %Compute ls and g.
        phasefilename = strcat(curfilename(1:end-12),'small',curfilename(end-8:end));
        texidxname = strcat(curfilename(1:end-12),'texidx',curfilename(end-8:end));
        phasemap = single(imread(phasefilename));
        texidx = imread(texidxname); %map of texton index
%         [gx,gy]=gradient(phasemap);
%         gmag2 = gx.^2+gy.^2;
%         h=fspecial('gaussian',[min(5*param.glswstd,45) min(5*param.glswstd,45)],param.glswstd);
%         avggmag2 = imfilter(gmag2,h,'same');%Avg gradient mag
%         %Next, compute the spatial variance of the phase
%         meanphase = imfilter(phasemap,h,'same');%First, mean phase value
%         diffphase2 = (phasemap-meanphase).^2;%Phase difference
%         phasevar = imfilter(diffphase2,h,'same');%Phase variance
%         ls = 1./phasevar;
%         g = avggmag2./phasevar.^2;
        
        
        lblim = imread(curfilename);
        glandim = (lblim==1);
        glandimnohole = imfill(glandim,'holes');
        lblimnohole = lblim|glandimnohole;%Get all the glands+lumen inside it
        res.roisize = sum(sum(lblimnohole~=0));
       % figure(1);
       % imagesc(lblim);
        
        %Now extracting a thin layer of stroma surrounding each gland and
        %compute the ls and g around it
       
        stromaim = (lblim==2);
        stromaimerode = imerode(stromaim,strel('disk',param.stromawidth));
        stromastrand = stromaim - stromaimerode;
        
        glandexpanded = imdilate(glandim,strel('disk',param.cleftingwidth));
        glandstrand = glandexpanded - glandim;
        stromastrand=stromastrand.*glandstrand;
        %stromastrand = bwareaopen(stromastrand,2000);%Eliminate too small area
        res.stromastrand = stromastrand;
        res.texidx = texidx;
        stromaidx = find(stromastrand==1);
%         res.mean_g = mean(g(stromaidx));
%         res.mean_ls = mean(ls(stromaidx));
        res.hist_tex_idx = hist(texidx(stromaidx),[1:50])/length(stromaidx);%This is the histogram of texton in a stroma area surrounding the gland
        res.stromaarea = length(stromaidx);%Save the areas of the areas where we evaluate the stroma regions
        
        %Extracting areas around the glands and look at the histogram
        glandstrand2=imdilate(glandim,strel('disk',param.basalwidth))-imerode(glandim,strel('disk',param.basalwidth));
        glandstrand2 = imdilate(stromaim,strel('disk',param.basalwidth)).*glandstrand2;
        basalidx = find(glandstrand2==1);
       % figure(2);
       % imagesc(glandstrand2);drawnow;
        res.basal_hist_tex_idx = hist(texidx(basalidx),[1:50])/length(basalidx);
        
        
end
