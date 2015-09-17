%This program extracts the histogram of textons in for the regions that may
%contain the basal cells
%Author: Tan H. Nguyen
%Last update: May 5th, 2015
%University of Illinois at Urbana-Chamapaign
    clc;
    %clear all;
    close all;
    datafolder = 'H:\TMA_cores_and_diagnosis\diagnosis_of_vicky\';
    addpath('..\support');
    g3folder =strcat(datafolder,'g3\');
    nmfolder = strcat(datafolder,'nm\');
    texton_hist_folder = 'H:\TMA_cores_and_diagnosis\texdir\'
    g3files = dir(strcat(g3folder,'*_lbl_*g3.tif'));
    nmfiles = dir(strcat(nmfolder,'*_lbl_*nm.tif'));
    g3names = cell(0,1);
    nmnames = cell(0,1);
    n3 = size(g3files,1); %Number of G3 samples
    nm = size(nmfiles,1); %Number of nm samples
    param.smallestlumenarea=5000;
    param.minglandarea = 5000;
    ntextons = 50;
    
    stromawidtharr = [150];
    re_obtain_histogram_data = 0;
    nstromawidth = length(stromawidtharr);
    normaltexthist = zeros(0,ntextons);
    g3texthist = zeros(1,ntextons);
    numberofhistogramperroi = 100000;% The number of histogram we take per ROI
    if (re_obtain_histogram_data)
        for stromawidthidx = 1:nstromawidth
            stromawidth = stromawidtharr(stromawidthidx);
            disp(['Processing at width = ' num2str(stromawidth)])
            param.glswstd = 25; %Ls and g windows size
            param.stromawidth = stromawidth;
            param.cleftingwidth =90;
            param.basalwidth=10;
            ntextons = 50;
            textonfeat = zeros(0,ntextons);
            iter = 0 ;
            curn3=0;
            curnnm = 0;
            numericclass = zeros(0,1);
            g3basaltextonhist=zeros(0,ntextons);
            nmbasaltextonhist=zeros(0,ntextons);
                      
            for idx=1:max(n3,nm)
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
                    if (~exist(filename2save,'file'))
                         res=extract_histogram_from_a_single_core_g3_vs_nm(curg3filename,texton_hist_filename,roi_cord_filename,param);
                                               save(filename2save,'res'); %Save the current response
                    else
                        temp=load(filename2save,'res');
                        res = temp.res;
                    end
                    curnhist=size(res.basalhistdata,1);
                    nhist = min(curnhist,numberofhistogramperroi);
                    sampleidx = randperm(curnhist);
                    g3basaltextonhist(end+1:end+nhist,:) = res.basalhistdata(sampleidx(1:nhist),:);
                    curn3 = curn3+1;
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
                    if (~exist(filename2save,'file'))
                        res=extract_histogram_from_a_single_core_g3_vs_nm(curnmfilename,texton_hist_filename,roi_cord_filename,param);
                        save(filename2save,'res'); %Save the current response
                    else
                        temp=load(filename2save,'res');
                        res = temp.res;
                    end
                    curnhist=size(res.basalhistdata,1);
                    nhist = min(curnhist,numberofhistogramperroi);
                    sampleidx = randperm(curnhist);                   
                    nmbasaltextonhist(end+1:end+nhist,:) = res.basalhistdata(sampleidx(1:nhist),:);
                    curnnm = curnnm+1;
                    tres = toc; disp(['Extraction time ' num2str(tres)]);

                end
            end
            save(strcat(texton_hist_folder,'basal_texton_hist_60_.mat'),'nmbasaltextonhist','g3basaltextonhist','-v7.3')
        end
    end
    
    %Next, load the data and perform LDA analysis
    res = load(strcat(texton_hist_folder,'basal_texton_hist_60_.mat'),'nmbasaltextonhist','g3basaltextonhist');
    nmdata = res.nmbasaltextonhist;
    g3data = res.g3basaltextonhist;
    nmnsamples = size(nmdata,1);
    g3nsamples = size(g3data,1);
    nmmu = mean(nmdata,1);nmmu = nmmu(:);
    g3mu = mean(g3data,1);g3mu = g3mu(:);
    mu = nmmu*nmnsamples + g3mu*g3nsamples;
    sb = nmnsamples*(nmmu-mu)*(nmmu-mu)'+g3nsamples*(g3mu-mu)*(g3mu-mu)'; %See Hoiem's slides
    nmcov = cov(nmdata);
    g3cov = cov(g3data);
    
    
    
    
    
