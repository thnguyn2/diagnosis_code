%function diagnosis_demo_based_on_core
%This function performs the analysis for the reactive stroma. It is
%mentioned that the normal tissue has more smooth muscle than the 
%Author: Tan H. Nguyen
%Last update: March 8th, 2015
%University of Illinois at Urbana-Chamapaign
    clc;
    %clear all;
    close all;
    datafolder = 'E:\TMA_cores_and_diagnosis\diagnosis_of_vicky\';
    addpath('..\support');
    g3folder =strcat(datafolder,'g3\');
    nmfolder = strcat(datafolder,'nm\');
    g3files = dir(strcat(g3folder,'*_lbl_*g3.tif'));
    nmfiles = dir(strcat(nmfolder,'*_lbl_*nm.tif'));
    g3names = cell(0,1);
    nmnames = cell(0,1);
    n3 = size(g3files,1); %Number of G3 samples
    nm = size(nmfiles,1); %Number of nm samples
    param.smallestlumenarea=5000;
    param.minglandarea = 5000;
    ntextons = 50;
    
    %stromawidtharr = [3 5 10 15 20 25 30 35 40 50 60 70 80 90 100 110 120 130];
    stromawidtharr = [150];
    obtain_auc_array = 1;
    nstromawidth = length(stromawidtharr);
    if (obtain_auc_array)
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
            g3textonfeat=zeros(0,ntextons);
            nmtextonfeat=zeros(0,ntextons);
            g3stromaarea = zeros(0,1);
            nmstromaarea = zeros(0,1);

            for idx=1:max(n3,nm)
           % for idx=43:43
             
                iter = iter+1;
                %Step 1: count the number of true glands and the number of lumen
                %each gland has
                %First, load the label map and cut out small fused regions. Make
                %sure that a slight fusion will not be detected as fused glands.
                %This seems to have been done in the code that generate the label
                %map.
                if (idx<=n3)
                    curg3filename = strcat(g3folder,g3files(idx).name);
                    res=process_single_core_g3_vs_nm(curg3filename,param);
                    disp(['G3 - ' g3files(idx).name]);
                    slashpos = strfind(g3files(idx).name,'_');
                    if ((sum(isnan(res.basal_hist_tex_idx))>0)|(res.stromaarea<20000)|(res.roisize<0.5*1e+6))
                        movefile(curg3filename,strcat(g3folder,'confused core\',g3files(idx).name));%Move the file to the new new location
                        %error(strcat(curg3filename,'has NaN histogram'));
                    else 
                        textonfeat(end+1,:)=res.basal_hist_tex_idx(:)'; %Just look at the histogram of texton indices
                        g3textonfeat(end+1,:)=res.basal_hist_tex_idx(:)';
                        curn3 = curn3+1;
                        numericclass(end+1,:)=3;

                    end
                    g3names{end+1,1}=g3files(idx).name;
                  
                end

                if (idx<=nm)
                    curnmfilename = strcat(nmfolder,nmfiles(idx).name);
                    res=process_single_core_g3_vs_nm(curnmfilename,param);
                     disp(['NM - ' nmfiles(idx).name]);
                    slashpos = strfind(nmfiles(idx).name,'_');
                   
                    if ((sum(isnan(res.basal_hist_tex_idx))>0)|(res.stromaarea<20000)|(res.roisize<0.5*1e+6))
                        %error(strcat(curnmfilename,'has NaN histogram'));
                        movefile(curnmfilename,strcat(nmfolder,'confused core\',nmfiles(idx).name));%Move the file to the new new location
                    else
                        textonfeat(end+1,:)=res.basal_hist_tex_idx(:)';
                        nmtextonfeat(end+1,:)=res.basal_hist_tex_idx(:)';
                        curnnm = curnnm+1;
                        numericclass(end+1,:)=0;
                    end
                    nmnames{end+1,1}=nmfiles(idx).name;
                end

                class=numericclass==0;%Convert into true and false datatype
                feat = textonfeat;

                if (mod(idx,4)==1)
                    %Next, try to spot out the texton that gives highest AUC for class
                    %3 and nm and determine where to look for them
                    texton_auc_arr = zeros(ntextons,1);
                    for textonidx=1:ntextons
                            curfeat = feat(:,end+textonidx-ntextons);
                            b1 = glmfit(curfeat,class,'normal');%Generalized linear model
                            p1 = glmval(b1,curfeat,'logit');%Compute the fitted probability
                            [x1,y1,t,texton_auc_arr(textonidx)] = perfcurve(class,p1,true);%Compute the ROC curve by providing the groundtruth
                    end

                    figure(7);
                    subplot(121);bar(texton_auc_arr);xlabel('Texton indices');ylabel('AUC values');drawnow;
                    %find 5 textons that has highest auc and display them
                    [maxval,tex_idx_sorted]=sort(texton_auc_arr,'descend');
                end
                %Compute the mean and the std of two classes
                mean_g3 = mean(g3textonfeat,1);
                std_g3 = std(g3textonfeat,0,1);
                mean_nm = mean(nmtextonfeat,1);
                std_nm = std(nmtextonfeat,0,1);
                figure(7);
                subplot(122);
                errorbar([1:ntextons],mean_g3,std_g3,'or');
                hold on;
                errorbar([1:ntextons],mean_nm,std_nm,'^b');
                hold off;drawnow;
                legend('g3','nm');
                title('Mean and standard deviation of the texton histogram on stroma regions')

%                 figure(4);
%                 plot(g3textonfeat(:,9),g3textonfeat(:,31),'^r');hold on;
%                 plot(nmtextonfeat(:,9),nmtextonfeat(:,31),'ob');hold off;
%                 xlabel('Texton 9');ylabel('Texton 31');
%                 
% 
%                 [g3stromaarea_sort,sortidx] = sort(g3stromaarea,'descend');
%                 g3stromaarea =  g3stromaarea_sort;  
%                 g3textonfeat = g3textonfeat(sortidx,:);
%                 %g3names = g3names{sortidx,1};
%                 [nmstromaarea_sort,sortidx] = sort(nmstromaarea,'descend');
%                 nmstromaarea = nmstromaarea_sort;
%                 nmtextonfeat = nmtextonfeat(sortidx,:);
%                 %nmnames = nmnames{sortidx,1};
%                 figure(3);
%                 imagesc(g3textonfeat,[0 0.4]);title('G3 carcinoma');colorbar
%                 set(gca,'YTickLabel',g3names);
%                 set(gca,'YTick',1:size(g3textonfeat,1));
%                 figure(6);
%                 imagesc(nmtextonfeat,[0 0.4]);title('Normal');colorbar;
%                 set(gca,'YTickLabel',nmnames);
%                 set(gca,'YTick',1:size(nmtextonfeat,1));
%                 
%                 figure(7);
%                 subplot(121);
%                 imagesc(g3stromaarea);title('G3 carcinoma')
%                 set(gca,'YTickLabel',g3names);
%                 set(gca,'YTick',1:size(g3textonfeat,1));
%                 subplot(122);
%                 imagesc(nmstromaarea);title('Normal')
%                 set(gca,'YTickLabel',nmnames);
%                 set(gca,'YTick',1:size(nmtextonfeat,1));
%                 drawnow;
%                 
%                 
%                 
%                 g3ratio = g3textonfeat(:,29)./g3textonfeat(:,31);
%                 nmratio = nmtextonfeat(:,29)./nmtextonfeat(:,31);
%                 x=linspace(0.5,5,20);
%                 g3hist = hist(g3ratio,x);
%                 nmhist = hist(nmratio,x);
%                 figure(8);
%                 stairs(x,g3hist,'b');hold on;
%                 stairs(x,nmhist,'r');hold on;
%                 
                
                drawnow;


            end
            clc;
            save(strcat('cancer_vs_normal_texton_feat_around_stroma_dist_',num2str(stromawidth),'.mat'),...
                'nmtextonfeat','g3textonfeat','texton_auc_arr','mean_g3','mean_nm','std_g3','std_nm','-v7.3');
        end
    end
    
    obtain_auc_array = 0;
    nstromawidth = length(stromawidtharr);
    
    %First, just load the auc for different cores and display the data
    auc_arr_all_scale = zeros(nstromawidth,ntextons);
    colorarr = 'rgbmk';
    patternarr = 'o^*v';
    figure(1);hold off;
    for stromaidx = 1:nstromawidth
        stromawidth =  stromawidtharr(stromaidx);
        load(strcat('cancer_vs_normal_texton_feat_around_stroma_dist_',num2str(stromawidth),'.mat'));
        auc_arr_all_scale(stromaidx,:)=texton_auc_arr;
        plot(texton_auc_arr,strcat(patternarr(floor(stromaidx/5)+1),colorarr(mod(stromaidx,5)+1)));
        hold on;
        
    end
    legend(num2str(stromawidtharr(:)));
    xlabel('Texton index');
    ylabel('AUC');
    
   first_text_idx=36;
    second_text_idx=50;
    %Draw the features for two most important features
    g3_feat_1 = g3textonfeat(:,first_text_idx);
    g3_feat_2 = g3textonfeat(:,second_text_idx);
    nm_feat_1 = nmtextonfeat(:,first_text_idx);
    nm_feat_2 = nmtextonfeat(:,second_text_idx);
    figure(2);
    plot(g3_feat_1,g3_feat_2,'^r');
    hold on;
    plot(nm_feat_1,nm_feat_2,'og');
    xlabel('Frequency of texton 1');
    ylabel('Frequency of texton 2');
    grid on;title('Scatter plot of histogram of two most important textons');
    legend('G3','NM')
    
    %Draw the ROC of the most important textons
    b1 = glmfit(feat(:,[first_text_idx second_text_idx]),class,'normal');%Generalized linear model
    p1 = glmval(b1,feat(:,[first_text_idx second_text_idx]),'logit');%Compute the fitted probability
    [x,y,t,auc] = perfcurve(class,p1,true);%Compute the ROC curve by providing the groundtruth
    figure(3);
    plot(x,y,'-*b');
    title('ROC curve for the 2 most important features');
    
          
    %Now, generate the plot for location of important textons
    for idx=1:max(n3,nm)
        if (idx<=n3)
            
                    curg3filename = strcat(g3folder,g3files(idx).name);
                    lblim = imread(curg3filename);
                    glandmask = lblim==1;
                    phasefilename = strcat(curg3filename(1:end-12),'small',curg3filename(end-8:end));
                    texidxname = strcat(curg3filename(1:end-12),'texidx',curg3filename(end-8:end));
                    phasemap = single(imresize(imread(phasefilename),[size(lblim,1) size(lblim,2)]));
                    texidx = imread(texidxname); %map of texton index
                    imptextidx1 = (texidx==first_text_idx).*glandmask;
                    imptextidx1 = imdilate(imptextidx1,strel('disk',2));
                    imptextidx2 = (texidx==second_text_idx).*glandmask;
                    imptextidx2 = imdilate(imptextidx2,strel('disk',2));
                    [nrows,ncols]=size(phasemap);
                    rgbmap = zeros(nrows,ncols,3);
                    rgbmap(:,:,1)=cast(imptextidx1,'single');
                    %rgbmap(:,:,2)=cast(imptextidx2,'single');
                    rgbmap(:,:,2)=phasemap/65535;
                    rgbmap(:,:,3)=phasemap/65535;
                    figure(4);
                    subplot(121);
                    imagesc(rgbmap);drawnow;
                    title('G3 Carcinoma');
                 

                end

                if (idx<=nm)
                    curnmfilename = strcat(nmfolder,nmfiles(idx).name);
                    lblim = imread(curnmfilename);
                    glandmask = lblim==1;
                    phasefilename = strcat(curnmfilename(1:end-12),'small',curnmfilename(end-8:end));
                    texidxname = strcat(curnmfilename(1:end-12),'texidx',curnmfilename(end-8:end));
                    phasemap = single(imresize(imread(phasefilename),[size(lblim,1) size(lblim,2)]));
                    texidx = imread(texidxname); %map of texton index
                    imptextidx1 = (texidx==first_text_idx).*glandmask;
                    imptextidx1 = imdilate(imptextidx1,strel('disk',2));
                    imptextidx2 = (texidx==second_text_idx).*glandmask;
                    imptextidx2 = imdilate(imptextidx2,strel('disk',2));
                    
                    [nrows,ncols]=size(phasemap);
                    rgbmap = zeros(nrows,ncols,3);
                    rgbmap(:,:,1)=cast(imptextidx1,'single');
                    %rgbmap(:,:,2)=cast(imptextidx2,'single');
                    rgbmap(:,:,2)=phasemap/65535;
                    rgbmap(:,:,3)=phasemap/65535;
                    figure(4);
                    subplot(122);
                    imagesc(rgbmap);drawnow;
                    title('Normal');
                    pause
                end
    end     
    
%end

