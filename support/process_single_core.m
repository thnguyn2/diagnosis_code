function res=process_single_core(curfilename,param)
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

        computelsandg = 1;
        displayinfooneachgland=1;
        if (computelsandg)
            %Compute ls and g.
            phasefilename = strcat(curfilename(1:end-12),'small',curfilename(end-8:end));
            phasemap = cast(imread(phasefilename),'single');
            phasemap = phasemap*4.0/65536-0.5;%Convert to phase
            [gx,gy]=gradient(phasemap); %Comute the gradient information
            gmag2 = (gx.^2+gy.^2)/(14)^2;%The factor in the back is 14 pixels/microns
            h=fspecial('gaussian',[min(5*param.glswstd,40) min(5*param.glswstd,40)],param.glswstd);
            h=h/sum(h);
            hzpad = zeros(size(phasemap));
            hzpad(1:size(h,1),1:size(h,2))=h;
            hzpad=circshift(hzpad,[-floor(size(h,1)/2) -floor(size(h,1)/2)]);
            hzpadf = fft2(hzpad);
            tic;
            avggmag2 = real(ifft2(fft2(gmag2).*hzpadf));
            meanphase = real(ifft2(fft2(phasemap).*hzpadf));
            diffphase2 = (phasemap-meanphase).^2;%Phase difference
            phasevar = real(ifft2(diffphase2.*hzpadf));
            ls = 1./phasevar;
            ko = 2*pi/0.57; %0.57 is the wavelength in microns
            g = 1-avggmag2./(2*ko^2*(phasevar.^2+1e-6));
            g = g.*(g>=0);
            
        end
        
       
          
        %Now, start computing the feature
        smallestlumenarea=param.smallestlumenarea;
        lblim = imread(curfilename);
        glandim = (lblim==1);
        glandimnohole = imfill(glandim,'holes');
        glandboundary = imdilate(glandimnohole,strel('disk',20))-glandimnohole;
        glandboundary = imresize(glandboundary,size(phasemap),'nearest');
        glandboundaryidx = find(glandboundary==1);
        
        lblimnohole = lblim|glandimnohole;
        res.roisize = sum(sum(lblimnohole~=0));
        
        
        %This images is for displaying the computing area for glands
        phase_rgb = zeros(size(phasemap,1),size(phasemap,2),3);
        phasemap1 = phasemap;
        phasemap1(glandboundaryidx)=max(phasemap(:));
        phase_rgb(:,:,1)=phasemap1;
        phasemap1(glandboundaryidx)=min(phasemap(:));
        phase_rgb(:,:,2)=phasemap1;
        phase_rgb(:,:,3)=phasemap1;
        
        
        glandareas = regionprops(glandim,'Area');
        centroids = regionprops(glandim,'Centroid'); %Get the coordinate of all centroids
        perimeters = regionprops(glandim,'Perimeter');
        solidity = regionprops(glandim,'Solidity');
        equidiameters = regionprops(glandim,'EquivDiameter');
        glandpixidxlist = regionprops(glandim,'PixelIdxList');
        nglands = size(glandareas,1);
        %Spot out the gland areas that are completely surrounded by gland
        disp(strcat('Found: ',num2str(nglands),' glands'));        
        %Calculate the lumen area of each glands
        glandimages = regionprops(glandim,'Image');
        nlumen_arr = zeros(nglands,1);
        solid_arr = zeros(nglands,1);
        cir_arr = zeros(nglands,1);
        area_arr = zeros(nglands,1); %This is an array for areas of the glands
        distortion_arr = zeros(nglands,1); %Distortion aray = perimeter of the gland/(pi*gland equivalent diameter)
        fractal_arr = zeros(nglands,1);
        
        for glandidx=1:nglands
            solid_arr(glandidx,1)=solidity(glandidx).Solidity;
            area_arr(glandidx,1)=glandareas(glandidx).Area;
            curglandimage = glandimages(glandidx,1).Image;
            distortion_arr(glandidx,1)=perimeters(glandidx).Perimeter/equidiameters(glandidx).EquivDiameter/pi;
            cir_arr(glandidx,1)=perimeters(glandidx).Perimeter^2/(4*pi*glandareas(glandidx).Area);
            %Fill out the lumen area
            glandfill = imfill(curglandimage,'holes');
            fractal_val=BoxCountfracDim(glandfill);
            lumenmap = im2bw(glandfill - curglandimage);
            nlumen = size(regionprops(lumenmap,'Area'),1);%Get the number of lumen regions inside the gland
            nlumen_arr(glandidx)=nlumen;
            fractal_arr(glandidx)=fractal_val;
            if (displayinfooneachgland)
                %Display the number of lumen for each image
                cencoord=centroids(glandidx).Centroid;
                xcoord = cencoord(1);
                ycoord = cencoord(2);
            end
        end        
        res.nlumen_arr=nlumen_arr;
        res.distortion_arr = distortion_arr;
        res.cir_arr = cir_arr;
        res.area_arr=area_arr;
        res.mean_area=mean(area_arr);
        res.std_area = std(area_arr);
        res.ratio_area=res.std_area/res.mean_area;
        res.nglands = nglands;
        %res.mean_nlum = sum(nlumen_arr.*area_arr)/sum(area_arr);
        res.mean_nlum = mean(nlumen_arr);
        res.max_nlum = max(nlumen_arr(:));
        res.mean_fractal = mean(fractal_arr);
        res.std_fractal = std(fractal_arr);
        %Compute the portion of fused glands
        total_fused_area = sum(area_arr(find(nlumen_arr>=2)));
        total_non_fused_area = sum(area_arr(find(nlumen_arr<=1)));
        res.fused_ratio = total_fused_area/(total_fused_area+total_non_fused_area);
      
        %res.fused_ratio = total_fused_area/sum(sum(lblim~=0));
      
        nlum_hist_bin = linspace(1,5,param.nbins);
        res.hist_nlum = hist(nlumen_arr,nlum_hist_bin);%Create a histogram of number of lumen
        res.hist_nlum = res.hist_nlum/sum(res.hist_nlum);
        res.hist_nlum = res.hist_nlum(:)';
        
        res.mean_dist = sum(distortion_arr.*area_arr)/sum(area_arr);
        distort_hist_bin = linspace(1,4,param.nbins);
        res.mean_dist = mean(distortion_arr);
        res.hist_dist = hist(distortion_arr,distort_hist_bin);
        res.hist_dist = res.hist_dist/sum(res.hist_dist);
        res.hist_dist = res.hist_dist(:)';
        res.med_dist = median(distortion_arr);
        res.mean_sold =  mean(solid_arr);
        
        res.mean_cir = mean(cir_arr);
        res.med_cir = median(cir_arr);   
        res.g_val_mean = mean(g(glandboundaryidx));
        res.g_val_median = median(g(glandboundaryidx));
        
        
        
        averageTime = toc/1000;
        disp(['Average extracting time' num2str(averageTime)]);

        
end
