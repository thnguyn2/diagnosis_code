function extract_histogram_of_texton_in_boundary_region
%This program extracts the histogram of textons in for the regions that may
%contain the basal cells
%Author: Tan H. Nguyen
%Last update: Oct 15th, 2015
%University of Illinois at Urbana-Chamapaign
    clc;
    clear all;
    %close all;
    datafolder = 'V:/TMA_cores_and_diagnosis/diagnosis_of_vicky/';
    addpath('../support');
    g3folder =strcat(datafolder,'g3/');
    g4folder =strcat(datafolder,'g4/');
    nmfolder = strcat(datafolder,'nm/');
    hgpfolder =strcat(datafolder,'hgp/');
    bphfolder = strcat(datafolder,'bph/');

    texton_hist_folder = 'V:/TMA_cores_and_diagnosis/texdir/'
    g3files = dir(strcat(g3folder,'*_lbl_*g3.tif'));
    g4files = dir(strcat(g4folder,'*_lbl_*g4.tif'));
    nmfiles = dir(strcat(nmfolder,'*_lbl_*nm.tif'));
    hgpfiles = dir(strcat(hgpfolder,'*_lbl_*hgp.tif'));
    bphfiles = dir(strcat(bphfolder,'*_lbl_*bph.tif'));
    g3names = cell(0,1);
    g4names = cell(0,1);
    nmnames = cell(0,1);
    hgpnames = cell(0,1);
    bphnames = cell(0,1);
    n3 = size(g3files,1); %Number of G3 samples
    n4 = size(g4files,1); %Number of G3 samples
    nm = size(nmfiles,1); %Number of nm samples
    nhgp = size(hgpfiles,1); %Number of G3 samples
    nbph = size(bphfiles,1); %Number of nm samples
    param.smallestlumenarea=5000;
    param.minglandarea = 5000;
    ntextons = 50;
    
    stromawidtharr = [200];
    re_obtain_histogram_data = 0;
    nstromawidth = length(stromawidtharr);
    normaltexthist = zeros(0,ntextons);
    g3texthist = zeros(0,ntextons);
    g4texthist = zeros(0,ntextons);
    bphtexthist = zeros(0,ntextons);
    hgptexthist = zeros(0,ntextons);
    
    numberofhistogramperroi = 80000;% The number of histogram we take per ROI
    reupdatebasalhist = 0;
    if (re_obtain_histogram_data)
        for stromawidthidx = 1:nstromawidth
            stromawidth = stromawidtharr(stromawidthidx);
            disp(['Processing at width = ' num2str(stromawidth)])
            param.glswstd = 25; %Ls and g windows size
            param.stromawidth = stromawidth;
            param.cleftingwidth =90;
            param.basalwidth=50;
            ntextons = 50;
            textonfeat = zeros(0,ntextons);
            iter = 0 ;
            curn4 = 0;
            curn3=0;
            curnnm = 0;
            curnbph=0;
            curnhgp = 0;
            numericclass = zeros(0,1);
            g4basaltextonhist=zeros(0,ntextons);
            g3basaltextonhist=zeros(0,ntextons);
            nmbasaltextonhist=zeros(0,ntextons);
            bphbasaltextonhist=zeros(0,ntextons);
            hgpbasaltextonhist=zeros(0,ntextons);
            g4glandtextonhist=zeros(0,ntextons);
            g3glandtextonhist=zeros(0,ntextons);
            nmglandtextonhist=zeros(0,ntextons);
            bphglandtextonhist=zeros(0,ntextons);
            hgpglandtextonhist=zeros(0,ntextons);
                      
            for idx=1:max([n3 n4 nm nbph nhgp])
                iter = iter+1;
                %Step 1: count the number of true glands and the number of lumen
                %each gland has
                %First, load the label map and cut out small fused regions. Make
                %sure that a slight fusion will not be detected as fused glands.
                %This seems to have been done in the code that generate the label
                %map.
                if (idx<=n3)
                    curg3filename = strcat(g3folder,g3files(idx).name);
                    dash_pos = find(g3files(idx).name=='_');
                    roi_cord_filename = strcat(g3folder,...
                        g3files(idx).name(1:dash_pos(1)),'coord',g3files(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI
                    roi_cord_filename(end-2:end)='mat';
                    curlabelname = g3files(idx).name(1:(dash_pos(1)-1)); %Label of the roi
                    texton_hist_filename=strcat(texton_hist_folder,curlabelname,'_texton_hist_60.mat');
                    tic;
                    filename2save = strcat(g3folder,...
                        g3files(idx).name(1:dash_pos(1)),'basal_texton_hist',g3files(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI)
                    filename2save(end-2:end)='mat';
                
                    disp(['G3 - ' g3files(idx).name]);
                    if ((~exist(filename2save,'file'))|( reupdatebasalhist))
                         res=extract_histogram_from_a_single_core_g3_vs_nm(curg3filename,texton_hist_filename,roi_cord_filename,param);
                         save(filename2save,'res','-v7.3'); %Save the current response
                    else
                        temp=load(filename2save,'res');
                        res = temp.res;
                    end
                    curnhist=size(res.basalhistdata,1);
                    nhist = min(curnhist,numberofhistogramperroi);
                    sampleidx = randperm(curnhist);
                    g3basaltextonhist(end+1:end+nhist,:) = res.basalhistdata(sampleidx(1:nhist),:);
                    
                    curnhist=size(res.glandhistdata,1);
                    nhist = min(curnhist,numberofhistogramperroi);
                    sampleidx = randperm(curnhist);
                    g3glandtextonhist(end+1:end+nhist,:) = res.glandhistdata(sampleidx(1:nhist),:);
                    
                    curn3 = curn3+1;
                    tres = toc; disp(['Extraction time ' num2str(tres)]);

                end
                
                if (idx<=n4)
                    curg4filename = strcat(g4folder,g4files(idx).name);
                    dash_pos = find(g4files(idx).name=='_');
                    roi_cord_filename = strcat(g4folder,...
                        g4files(idx).name(1:dash_pos(1)),'coord',g4files(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI
                    roi_cord_filename(end-2:end)='mat';
                    curlabelname = g4files(idx).name(1:(dash_pos(1)-1)); %Label of the roi
                    texton_hist_filename=strcat(texton_hist_folder,curlabelname,'_texton_hist_60.mat');
                    tic;
                    filename2save = strcat(g4folder,...
                        g4files(idx).name(1:dash_pos(1)),'basal_texton_hist',g4files(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI)
                    filename2save(end-2:end)='mat';
                
                    disp(['G4 - ' g4files(idx).name]);
                    if ((~exist(filename2save,'file'))|(reupdatebasalhist))
                         res=extract_histogram_from_a_single_core_g3_vs_nm(curg4filename,texton_hist_filename,roi_cord_filename,param);
                         save(filename2save,'res','-v7.3'); %Save the current response
                    else
                        temp=load(filename2save,'res');
                        res = temp.res;
                    end
                    curnhist=size(res.basalhistdata,1);
                    nhist = min(curnhist,numberofhistogramperroi);
                    sampleidx = randperm(curnhist);
                    g4basaltextonhist(end+1:end+nhist,:) = res.basalhistdata(sampleidx(1:nhist),:);
                    
                    curnhist=size(res.glandhistdata,1);
                    nhist = min(curnhist,numberofhistogramperroi);
                    sampleidx = randperm(curnhist);
                    g4glandtextonhist(end+1:end+nhist,:) = res.glandhistdata(sampleidx(1:nhist),:);
                    
                    curn4 = curn4+1;
                    tres = toc; disp(['Extraction time ' num2str(tres)]);

                end

                if (idx<=nm)
                    curnmfilename = strcat(nmfolder,nmfiles(idx).name);
                    dash_pos = find(nmfiles(idx).name=='_');
                    roi_cord_filename = strcat(nmfolder,...
                    nmfiles(idx).name(1:dash_pos(1)),'coord',nmfiles(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI
                    roi_cord_filename(end-2:end)='mat';
                    curlabelname = nmfiles(idx).name(1:(dash_pos(1)-1)); %Label of the roi
                    texton_hist_filename=strcat(texton_hist_folder,curlabelname,'_texton_hist_60.mat');
                    tic;
                    filename2save = strcat(nmfolder,...
                        nmfiles(idx).name(1:dash_pos(1)),'basal_texton_hist',nmfiles(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI)
                    filename2save(end-2:end)='mat';
                    disp(['NM - ' nmfiles(idx).name]);
                    if ((~exist(filename2save,'file'))|(reupdatebasalhist))
                        res=extract_histogram_from_a_single_core_g3_vs_nm(curnmfilename,texton_hist_filename,roi_cord_filename,param);
                        save(filename2save,'res','-v7.3'); %Save the current response
                    else
                        temp=load(filename2save,'res');
                        res = temp.res;
                    end
                    curnhist=size(res.basalhistdata,1);
                    nhist = min(curnhist,numberofhistogramperroi);
                    sampleidx = randperm(curnhist);                   
                    nmbasaltextonhist(end+1:end+nhist,:) = res.basalhistdata(sampleidx(1:nhist),:);
                    
                    curnhist=size(res.glandhistdata,1);
                    nhist = min(curnhist,numberofhistogramperroi);
                    sampleidx = randperm(curnhist);                   
                    nmglandtextonhist(end+1:end+nhist,:) = res.glandhistdata(sampleidx(1:nhist),:);
                    
                    curnnm = curnnm+1;
                    tres = toc; disp(['Extraction time ' num2str(tres)]);

                end
                
                 if (idx<=nbph)
                    curbphfilename = strcat(bphfolder,bphfiles(idx).name);
                    dash_pos = find(bphfiles(idx).name=='_');
                    roi_cord_filename = strcat(bphfolder,...
                        bphfiles(idx).name(1:dash_pos(1)),'coord',bphfiles(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI
                    roi_cord_filename(end-2:end)='mat';
                    curlabelname = bphfiles(idx).name(1:(dash_pos(1)-1)); %Label of the roi
                    texton_hist_filename=strcat(texton_hist_folder,curlabelname,'_texton_hist_60.mat');
                    tic;
                    filename2save = strcat(bphfolder,...
                        bphfiles(idx).name(1:dash_pos(1)),'basal_texton_hist',bphfiles(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI)
                    filename2save(end-2:end)='mat';
                
                    disp(['BPH - ' bphfiles(idx).name]);
                    if ((~exist(filename2save,'file'))|(reupdatebasalhist))
                         res=extract_histogram_from_a_single_core_g3_vs_nm(curbphfilename,texton_hist_filename,roi_cord_filename,param);
                         save(filename2save,'res','-v7.3'); %Save the current response
                    else
                        temp=load(filename2save,'res');
                        res = temp.res;
                    end
                    curnhist=size(res.basalhistdata,1);
                    nhist = min(curnhist,numberofhistogramperroi);
                    sampleidx = randperm(curnhist);
                    bphbasaltextonhist(end+1:end+nhist,:) = res.basalhistdata(sampleidx(1:nhist),:);
                    
                    curnhist=size(res.glandhistdata,1);
                    nhist = min(curnhist,numberofhistogramperroi);
                    sampleidx = randperm(curnhist);
                    bphglandtextonhist(end+1:end+nhist,:) = res.glandhistdata(sampleidx(1:nhist),:);
                    
                    curnbph = curnbph+1;
                    tres = toc; disp(['Extraction time ' num2str(tres)]);

                 end
                
                  if (idx<=nhgp)
                    curhgpfilename = strcat(hgpfolder,hgpfiles(idx).name);
                    dash_pos = find(hgpfiles(idx).name=='_');
                    roi_cord_filename = strcat(hgpfolder,...
                        hgpfiles(idx).name(1:dash_pos(1)),'coord',hgpfiles(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI
                    roi_cord_filename(end-2:end)='mat';
                    curlabelname = hgpfiles(idx).name(1:(dash_pos(1)-1)); %Label of the roi
                    texton_hist_filename=strcat(texton_hist_folder,curlabelname,'_texton_hist_60.mat');
                    tic;
                    filename2save = strcat(hgpfolder,...
                        hgpfiles(idx).name(1:dash_pos(1)),'basal_texton_hist',hgpfiles(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI)
                    filename2save(end-2:end)='mat';
                
                    disp(['hgp - ' hgpfiles(idx).name]);
                    if ((~exist(filename2save,'file'))|(reupdatebasalhist))
                         res=extract_histogram_from_a_single_core_g3_vs_nm(curhgpfilename,texton_hist_filename,roi_cord_filename,param);
                         save(filename2save,'res','-v7.3'); %Save the current response
                    else
                        temp=load(filename2save,'res');
                        res = temp.res;
                    end
                    curnhist=size(res.basalhistdata,1);
                    nhist = min(curnhist,numberofhistogramperroi);
                    sampleidx = randperm(curnhist);
                    hgpbasaltextonhist(end+1:end+nhist,:) = res.basalhistdata(sampleidx(1:nhist),:);
                    
                    curnhist=size(res.glandhistdata,1);
                    nhist = min(curnhist,numberofhistogramperroi);
                    sampleidx = randperm(curnhist);
                    hgpglandtextonhist(end+1:end+nhist,:) = res.glandhistdata(sampleidx(1:nhist),:);
                    
                    curnhgp = curnhgp+1;
                    tres = toc; disp(['Extraction time ' num2str(tres)]);

                end

            end
             save(strcat(texton_hist_folder,'basal_texton_hist_60_.mat'),'nmbasaltextonhist',...
                 'hgpbasaltextonhist','bphbasaltextonhist','g3basaltextonhist','g4basaltextonhist','nmglandtextonhist',...
                 'hgpglandtextonhist','bphglandtextonhist','g3glandtextonhist','g4glandtextonhist','-v7.3')
        end
    end
    
    recomputecov =0;
    nhistperclass = 20000000;
    benigndata = zeros(0,50);
    cancerdata = zeros(0,50);
    if ((~exist(strcat(texton_hist_folder,'basal_lda.mat')))|(recomputecov))
         %Next, load the data and perform LDA analysis
        res = load(strcat(texton_hist_folder,'basal_texton_hist_60_.mat'),'nmbasaltextonhist');        nmdata = res.nmbasaltextonhist;
        nmdata = nmdata(1:min(size(nmdata,1),nhistperclass),:);
        clear res;
        nmmu = mean(nmdata,1);nmmu = nmmu(:);
        nmcov = cov(nmdata);
        nmnsamples = size(nmdata,1);
        benigndata(end+1:end+nmnsamples,:) = nmdata;
        clear nmdata;
        
        
        res = load(strcat(texton_hist_folder,'basal_texton_hist_60_.mat'),'g3basaltextonhist');
        g3data = res.g3basaltextonhist;
        g3data = g3data(1:min(size(g3data,1),nhistperclass),:);
        clear res;
        g3mu = mean(g3data,1);g3mu = g3mu(:);
        g3cov = cov(g3data);
        g3nsamples = size(g3data,1);
        cancerdata(end+1:end+g3nsamples,:) = g3data;
        clear g3data;
        
        res = load(strcat(texton_hist_folder,'basal_texton_hist_60_.mat'),'g4basaltextonhist');
        g4data = res.g4basaltextonhist;
        g4data = g4data(1:min(size(g4data,1),nhistperclass),:);
        clear res;
        g4mu = mean(g4data,1);g3mu = g4mu(:);
        g4cov = cov(g4data);
        g4nsamples = size(g4data,1);
        cancerdata(end+1:end+g4nsamples,:) = g4data;
        clear g4data;
        
        res = load(strcat(texton_hist_folder,'basal_texton_hist_60_.mat'),'bphbasaltextonhist');
        bphdata = res.bphbasaltextonhist;
        bphdata = bphdata(1:min(size(bphdata,1),nhistperclass),:);
        clear res;
        bphmu = mean(bphdata,1);bphmu = bphmu(:);
        bphcov = cov(bphdata);
        bphnsamples = size(bphdata,1);
        benigndata(end+1:end+bphnsamples,:) = bphdata;
        clear bphdata;
        
        res = load(strcat(texton_hist_folder,'basal_texton_hist_60_.mat'),'hgpbasaltextonhist');
        hgpdata = res.hgpbasaltextonhist;
        hgpdata = hgpdata(1:min(size(hgpdata,1),nhistperclass),:);
        clear res;
        hgpmu = mean(hgpdata,1);hgpmu = hgpmu(:);
        hgpcov = cov(hgpdata);
        hgpnsamples = size(hgpdata,1);
        benigndata(end+1:end+hgpnsamples,:) = hgpdata;
        clear hgpdata;
        
        benignmu = mean(benigndata,1);
        cancermu = mean(cancerdata,1);
        benigncov = cov(benigndata);
        cancercov = cov(cancerdata);
        
        save(strcat(texton_hist_folder,'basal_lda.mat'),'benignmu','cancermu',...
            'benigncov','cancercov',...
            'g3mu','g4mu','nmmu','bphmu','hgpmu','g3cov','g4cov','nmcov','hgpcov','bphcov','bphnsamples','hgpnsamples','g4nsamples','g3nsamples','nmnsamples');
    else        load(strcat(texton_hist_folder,'basal_lda.mat'));
    end
    
%   sw = g3cov + benigncov;
%   sb = (g3mu-benignmu(:))*(g3mu-benignmu(:))';
    
     sw = benigncov + cancercov;
     sb = (benignmu(:)-cancermu(:))*(benignmu(:)-cancermu(:))';
    [v,d]=eig(inv(sw)*sb);
    
     recomputecov = 0;
     gbenigndata = zeros(0,50);
     gcancerdata = zeros(0,50);
     if ((~exist(strcat(texton_hist_folder,'gland_lda.mat')))|(recomputecov))%Compute the mean and cov matrix using the gland features
         %Next, load the data and perform LDA analysis
        res = load(strcat(texton_hist_folder,'basal_texton_hist_60_.mat'),'nmglandtextonhist');
        gnmdata = res.nmglandtextonhist;
        gnmdata = gnmdata(1:min(size(gnmdata,1),nhistperclass),:);
        clear res;
        gnmmu = mean(gnmdata,1);gnmmu = gnmmu(:);
        gnmcov = cov(gnmdata);
        gnmnsamples = size(gnmdata,1);
        gbenigndata(end+1:end+gnmnsamples,:) = gnmdata;
        clear gnmdata;

        res = load(strcat(texton_hist_folder,'basal_texton_hist_60_.mat'),'g3glandtextonhist');
        gg3data = res.g3glandtextonhist;
        gg3data = gg3data(1:min(size(gg3data,1),nhistperclass),:);
        clear res;
        gg3mu = mean(gg3data,1);gg3mu = gg3mu(:);
        gg3cov = cov(gg3data);
        gg3nsamples = size(gg3data,1);
        gcancerdata(end+1:end+gg3nsamples,:) = gg3data;
        clear gg3data;
      
        
        res = load(strcat(texton_hist_folder,'basal_texton_hist_60_.mat'),'g4glandtextonhist');
        gg4data = res.g4glandtextonhist;
        gg4data = gg4data(1:min(size(gg4data,1),nhistperclass),:);
        clear res;
        gg4mu = mean(gg4data,1);gg3mu = gg3mu(:);
        gg4cov = cov(gg4data);
        gg4nsamples = size(gg4data,1);
        gcancerdata(end+1:end+gg4nsamples,:) = gg4data;
        clear gg4data;

        res = load(strcat(texton_hist_folder,'basal_texton_hist_60_.mat'),'bphglandtextonhist');
        gbphdata = res.bphglandtextonhist;
        gbphdata = gbphdata(1:min(size(gbphdata,1),nhistperclass),:);
        clear res;
        gbphmu = mean(gbphdata,1);gbphmu = gbphmu(:);
        gbphcov = cov(gbphdata);
        gbphnsamples = size(gbphdata,1);
        gbenigndata(end+1:end+gbphnsamples,:) = gbphdata;
        clear gbphdata;
        
        res = load(strcat(texton_hist_folder,'basal_texton_hist_60_.mat'),'hgpglandtextonhist');
        ghgpdata = res.hgpglandtextonhist;
        ghgpdata = ghgpdata(1:min(size(ghgpdata,1),nhistperclass),:);
        clear res;
        ghgpmu = mean(ghgpdata,1);ghgpmu = ghgpmu(:);
        ghgpcov = cov(ghgpdata);
        ghgpnsamples = size(ghgpdata,1);
        gbenigndata(end+1:end+ghgpnsamples,:) = ghgpdata;
        clear ghgpdata;
        
        gbenignmu = mean(gbenigndata,1);
        gcancermu = mean(gcancerdata,1);
        gbenigncov = cov(gbenigndata);
        gcancercov = cov(gcancerdata);
        
        gsw = gg3cov + gnmcov;
        gsb = (gg3mu-gnmmu)*(gg3mu-gnmmu)';
        [gv,gd]=eig(inv(gsw)*gsb);
        save(strcat(texton_hist_folder,'gland_lda.mat'),'gbenignmu','gcancermu',...
            'gbenigncov','gcancercov','gv','gd','gg3mu','gg4mu','gnmmu','gbphmu','ghgpmu','gg4cov','gg3cov','gnmcov','ghgpcov','gbphcov',...
            'gbphnsamples','ghgpnsamples','gg4nsamples','gg3nsamples','gnmnsamples');
    else
        load(strcat(texton_hist_folder,'gland_lda.mat'));
    end
    
    
    gsw = gbenigncov + gcancercov;
    gsb = (gbenignmu(:)-gcancermu(:))*(gbenignmu(:)-gcancermu(:))';
    
    %gsw = gg3cov + gbenigncov;
    %gsb = (gg3mu-gbenignmu(:))*(gg3mu-gbenignmu(:))';
    [gv,gd]=eig(inv(gsw)*gsb);
 
     %Now, go over all the cores and draw the projection of the histogram.
     %The numeric value will show locations that is significantly different
     w=v(:,1);
     gw = gv(:,1);
     maxproj = 300;
     minproj = -300;
     for idx=1:max([n3 n4 nm nbph nhgp])
                if (idx<=n3)
                    curg3filename = strcat(g3folder,g3files(idx).name);
                    dash_pos = find(g3files(idx).name=='_');
                    filename2read = strcat(g3folder,...
                        g3files(idx).name(1:dash_pos(1)),'basal_texton_hist',g3files(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI)
                    filename2read(end-2:end)='mat';
                    phaseimname = strcat(g3folder,...
                        g3files(idx).name(1:dash_pos(1)),'small',g3files(idx).name(dash_pos(2):end));
                    
                    disp(['G3 - ' g3files(idx).name]);
                    temp=load(filename2read,'res');
                   if (isfield(temp,'res')) %Need to figure why it does not work for T14 G4

                        res = temp.res;
                        projdata = w'*res.basalhistdata';
                        gprojdata = w'*res.glandhistdata';                    
                        %positivepart = length(find(projdata>=10))/length(projdata);
                        %negativepart = length(find(projdata<=-10))/length(projdata);
                        %figure(2);hold on; plot(positivepart,negativepart,'*r');axis([0 1 0 1]);
                        figure(2);plot(mean(projdata),mean(gprojdata),'*r');hold on;
                        phasemap = imresize(single(imread(phaseimname)),[res.nrows,res.ncols],'nearest')/65536.0;
                        newmap=zeros(size(phasemap));
                        newmap(res.basalidx)=(projdata);
                        newmap(res.glandidx)=gprojdata;
                        H = convert2hsv(phasemap,newmap,maxproj,minproj); 
                        figure(1);
                        subplot(221);imshow(H);title('G3');drawnow;
                        filename2save = strcat(g3folder,...
                        g3files(idx).name(1:dash_pos(1)),'basal_heat_map',g3files(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI)
                        %writeTIFF(H,filename2save,'uint16');
                       
                   end
                end
                
                if (idx<=n4)
                    curg4filename = strcat(g4folder,g4files(idx).name);
                    dash_pos = find(g4files(idx).name=='_');
                    filename2read = strcat(g4folder,...
                        g4files(idx).name(1:dash_pos(1)),'basal_texton_hist',g4files(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI)
                    filename2read(end-2:end)='mat';
                    phaseimname = strcat(g4folder,...
                        g4files(idx).name(1:dash_pos(1)),'small',g4files(idx).name(dash_pos(2):end));
                    
                    disp(['G4 - ' g4files(idx).name]);
                    temp=load(filename2read,'res');
                    if (isfield(temp,'res')) %Need to figure why it does not work for T14 G4
                        res = temp.res;
                        projdata = w'*res.basalhistdata';
                        gprojdata = w'*res.glandhistdata';                    

                        figure(2);plot(mean(projdata),mean(gprojdata),'*k');hold on;
                        phasemap = imresize(single(imread(phaseimname)),[res.nrows,res.ncols])/65535.0;
                        newmap=zeros(size(phasemap));
                        newmap(res.basalidx)=(projdata);
                        newmap(res.glandidx)=gprojdata;

                        H = convert2hsv(phasemap,newmap,maxproj,minproj); 
                        figure(1);
                        subplot(221);imshow(H);title('G4');drawnow;
                        filename2save = strcat(g4folder,...
                            g4files(idx).name(1:dash_pos(1)),'basal_heat_map',g4files(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI)
                        %writeTIFF(H,filename2save,'uint16');
                       
                    end
                end

                if (idx<=nm)
                    curnmfilename = strcat(nmfolder,nmfiles(idx).name);
                    dash_pos = find(nmfiles(idx).name=='_');                    
                    filename2read = strcat(nmfolder,...
                        nmfiles(idx).name(1:dash_pos(1)),'basal_texton_hist',nmfiles(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI)
                    filename2read(end-2:end)='mat';
                    phaseimname = strcat(nmfolder,...
                        nmfiles(idx).name(1:dash_pos(1)),'small',nmfiles(idx).name(dash_pos(2):end));
                
                    disp(['NM - ' nmfiles(idx).name]);
                    if (exist(filename2read))
                        temp=load(filename2read,'res');
                        res = temp.res;
                        projdata = w'*res.basalhistdata';
                        gprojdata = w'*res.glandhistdata';                    
                 
                        figure(2);plot(mean(projdata),mean(gprojdata),'ob');hold on;
                        phasemap = imresize(single(imread(phaseimname)),[res.nrows,res.ncols])/65535.0;
                        newmap=zeros(size(phasemap));
                        newmap(res.basalidx)=(projdata);
                        newmap(res.glandidx)=gprojdata;
                        H = convert2hsv(phasemap,newmap,maxproj,minproj); 
                        figure(1);
                        subplot(222);imshow(H);title('Normal')

                        filename2save = strcat(nmfolder,...
                            nmfiles(idx).name(1:dash_pos(1)),'basal_heat_map',nmfiles(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI)
                        %writeTIFF(H,filename2save,'uint16'); 
                      
                    end
                end
                
                if (idx<=nhgp)
                    curhgpfilename = strcat(hgpfolder,hgpfiles(idx).name);
                    dash_pos = find(hgpfiles(idx).name=='_');
                    filename2read = strcat(hgpfolder,...
                        hgpfiles(idx).name(1:dash_pos(1)),'basal_texton_hist',hgpfiles(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI)
                    filename2read(end-2:end)='mat';
                    phaseimname = strcat(hgpfolder,...
                        hgpfiles(idx).name(1:dash_pos(1)),'small',hgpfiles(idx).name(dash_pos(2):end));
                    
                    disp(['HGP - ' hgpfiles(idx).name]);
                    temp=load(filename2read,'res');
                    res = temp.res;
                    projdata = w'*res.basalhistdata';
                    gprojdata = w'*res.glandhistdata';                    
                 
%                     positivepart = length(find(projdata>=100))/length(projdata);
%                     negativepart = length(find(projdata<=-10))/length(projdata);
%                     figure(2);hold on; plot(positivepart,negativepart,'og');axis([0 1 0 1]);
                    figure(2);plot(mean(projdata),mean(gprojdata),'og');hold on;
                    phasemap = imresize(single(imread(phaseimname)),[res.nrows,res.ncols])/65535.0;
                    newmap=zeros(size(phasemap));
                    newmap(res.basalidx)=(projdata);
                    newmap(res.glandidx)=gprojdata;
                    
                    H = convert2hsv(phasemap,newmap,maxproj,minproj); 
                    figure(1);
                    subplot(223);imshow(H);title('HGP');drawnow;
                    filename2save = strcat(hgpfolder,...
                        hgpfiles(idx).name(1:dash_pos(1)),'basal_heat_map',hgpfiles(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI
                    %writeTIFF(H,filename2save,'uint16'); 
                   
                end
                
                 if (idx<=nbph)
                    curbphfilename = strcat(bphfolder,bphfiles(idx).name);
                    dash_pos = find(bphfiles(idx).name=='_');                    
                    filename2read = strcat(bphfolder,...
                        bphfiles(idx).name(1:dash_pos(1)),'basal_texton_hist',bphfiles(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI)
                    filename2read(end-2:end)='mat';
                    phaseimname = strcat(bphfolder,...
                        bphfiles(idx).name(1:dash_pos(1)),'small',bphfiles(idx).name(dash_pos(2):end));
                
                    disp(['BPH - ' bphfiles(idx).name]);
                    if (exist(filename2read))
                        temp=load(filename2read,'res');
                        res = temp.res;
                        projdata = w'*res.basalhistdata';
                        gprojdata = w'*res.glandhistdata';                    
                 
%                         positivepart = length(find(projdata>=100))/length(projdata);
%                         negativepart = length(find(projdata<=-10))/length(projdata);
%                         figure(2);hold on; plot(positivepart,negativepart,'om');axis([0 1 0 1]);
                        figure(2);plot(mean(projdata),mean(gprojdata),'om');hold on;
                         phasemap = imresize(single(imread(phaseimname)),[res.nrows,res.ncols])/65535.0;
                       newmap=zeros(size(phasemap));
                        newmap(res.basalidx)=(projdata);
                        newmap(res.glandidx)=gprojdata;
                        
                        H = convert2hsv(phasemap,newmap,maxproj,minproj); 
                        figure(1);
                        subplot(224);imshow(H);title('BPH');drawnow;
                        filename2save = strcat(bphfolder,...
                            bphfiles(idx).name(1:dash_pos(1)),'basal_heat_map',nmfiles(idx).name(dash_pos(2):end));%Coordinates of the region to extract the ROI)
                        %writeTIFF(H,filename2save,'uint16'); 
                      
                    end
                end

     end

end


    
    
