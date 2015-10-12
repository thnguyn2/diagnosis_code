%This function performs the analysis for the tissue expression around each
%glands.It compute the histogram of texton indices around the glands at
%different window sizes.
%Author: Tan H. Nguyen The output is the histogram of texton data for all
%the cores
%Last update: Oct 11th, 2015
%University of Illinois at Urbana-Chamapaign
    clc;
    %clear all;
    close all;
    datafolder = 'V:/TMA_cores_and_diagnosis/diagnosis_of_vicky/';
    addpath('../support');
    g3folder =strcat(datafolder,'g3/');
    nmfolder = strcat(datafolder,'nm/');
    g4folder =strcat(datafolder,'g4/');
    bphfolder = strcat(datafolder,'bph/');
    hgpfolder = strcat(datafolder,'hgp/');
  
    param.smallestlumenarea=5000;
    param.minglandarea = 5000;
    ntextons = 50;
    
    stromawidtharr = [20 30 40 50 60 70 80 90 100 110 120];
    obtain_auc_array = 0;
    if (obtain_auc_array)
        g3files = dir(strcat(g3folder,'*_lbl_*g3.tif'));
        g4files = dir(strcat(g4folder,'*_lbl_*g4.tif'));
        nmfiles = dir(strcat(nmfolder,'*_lbl_*nm.tif'));
        hgpfiles = dir(strcat(hgpfolder,'*_lbl_*hgp.tif'));
        bphfiles = dir(strcat(bphfolder,'*_lbl_*bph.tif'));

        g3names = cell(0,1);
        nmnames = cell(0,1);
        g4names = cell(0,1);
        bphnames = cell(0,1);
        hgpnames = cell(0,1);

        n3 = size(g3files,1); %Number of G3 samples
        nm = size(nmfiles,1); %Number of nm samples
        n4 = size(g4files,1); %Number of G3 samples
        nbph = size(bphfiles,1); %Number of nm samples
        nhgp = size(hgpfiles,1); %Number of nm samples
        
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
            curn4 = 0;
            curnbph = 0;
            curnhgp = 0;
            numericclass = zeros(0,1);
            g3textonfeat=zeros(0,ntextons);
            g4textonfeat=zeros(0,ntextons);
            nmtextonfeat=zeros(0,ntextons);
            bphtextonfeat=zeros(0,ntextons);
            hgptextonfeat = zeros(0,ntextons);
            
            g3stromaarea = zeros(0,1);
            g4stromaarea = zeros(0,1);
            nmstromaarea = zeros(0,1);
            bphtromaarea = zeros(0,1);
            hgptromaarea = zeros(0,1);
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
                    res=process_single_core_g3_vs_nm(curg3filename,param,'g3');
                    disp(['G3 - ' g3files(idx).name]);
                    slashpos = strfind(g3files(idx).name,'_');
                    if ((sum(isnan(res.basal_hist_tex_idx))>0)|(res.stromaarea<20000)|(res.roisize<0.5*1e+6))
                        %movefile(curg3filename,strcat(g3folder,'confused core\',g3files(idx).name));%Move the file to the new new location
                    else 
                        textonfeat(end+1,:)=res.basal_hist_tex_idx(:)'; %Just look at the histogram of texton indices
                        g3textonfeat(end+1,:)=res.basal_hist_tex_idx(:)';
                        curn3 = curn3+1;
                        numericclass(end+1,:)=3;
                        g3names{end+1,1}=g3files(idx).name;
                    end
                    
                  
                end
                
                 if (idx<=n4)
                    curg4filename = strcat(g4folder,g4files(idx).name);
                    res=process_single_core_g3_vs_nm(curg4filename,param,'g4');
                    disp(['G4 - ' g4files(idx).name]);
                    slashpos = strfind(g4files(idx).name,'_');
                    if ((sum(isnan(res.basal_hist_tex_idx))>0)|(res.stromaarea<20000)|(res.roisize<0.5*1e+6))
                        %movefile(curg4filename,strcat(g4folder,'confused core\',g4files(idx).name));%Move the file to the new new location
                    else 
                        textonfeat(end+1,:)=res.basal_hist_tex_idx(:)'; %Just look at the histogram of texton indices
                        g4textonfeat(end+1,:)=res.basal_hist_tex_idx(:)';
                        curn4 = curn4+1;
                        numericclass(end+1,:)=4;
                        g4names{end+1,1}=g4files(idx).name;
                  
            
                    end
                end

                if (idx<=nm)
                    curnmfilename = strcat(nmfolder,nmfiles(idx).name);
                    if (exist(curnmfilename,'file'))
                        res=process_single_core_g3_vs_nm(curnmfilename,param,'nm');
                         disp(['NM - ' nmfiles(idx).name]);
                        slashpos = strfind(nmfiles(idx).name,'_');

                        if ((sum(isnan(res.basal_hist_tex_idx))>0)|(res.stromaarea<20000)|(res.roisize<0.5*1e+6))
                           %movefile(curnmfilename,strcat(nmfolder,'confused core\',nmfiles(idx).name));%Move the file to the new new location
                        else
                            textonfeat(end+1,:)=res.basal_hist_tex_idx(:)';
                            nmtextonfeat(end+1,:)=res.basal_hist_tex_idx(:)';
                            curnnm = curnnm+1;
                            numericclass(end+1,:)=0;
                            nmnames{end+1,1}=nmfiles(idx).name;

                        end
                    end
                end
                
                 if (idx<=nbph)
                    curbphfilename = strcat(bphfolder,bphfiles(idx).name);
                    if (exist(curbphfilename,'file'))
                        res=process_single_core_g3_vs_nm(curbphfilename,param,'bph');
                         disp(['BPH - ' bphfiles(idx).name]);
                        slashpos = strfind(bphfiles(idx).name,'_');

                        if ((sum(isnan(res.basal_hist_tex_idx))>0)|(res.stromaarea<20000)|(res.roisize<0.5*1e+6))
                           %movefile(curbphfilename,strcat(bphfolder,'confused core\',bphfiles(idx).name));%Move the file to the new new location
                        else
                            textonfeat(end+1,:)=res.basal_hist_tex_idx(:)';
                            nmtextonfeat(end+1,:)=res.basal_hist_tex_idx(:)';
                            curnbph = curnbph+1;
                            numericclass(end+1,:)=-1;
                            bphnames{end+1,1}=bphfiles(idx).name;

                        end
                    end
                end
                
                   if (idx<=nhgp)
                    curhgpfilename = strcat(hgpfolder,hgpfiles(idx).name);
                    if (exist(curhgpfilename,'file'))                      
                        res=process_single_core_g3_vs_nm(curhgpfilename,param,'hgp');
                         disp(['HGP - ' hgpfiles(idx).name]);
                        slashpos = strfind(hgpfiles(idx).name,'_');

                        if ((sum(isnan(res.basal_hist_tex_idx))>0)|(res.stromaarea<20000)|(res.roisize<0.5*1e+6))
                           %movefile(curhgpfilename,strcat(hgpfolder,'confused core\',hgpfiles(idx).name));%Move the file to the new new location
                        else
                            textonfeat(end+1,:)=res.basal_hist_tex_idx(:)';
                            hgptextonfeat(end+1,:)=res.basal_hist_tex_idx(:)';
                            curnhgp = curnhgp+1;
                            numericclass(end+1,:)=-2;
                            hgpnames{end+1,1}=hgpfiles(idx).name;

                        end
                    end
               end
                class=numericclass==0;%Convert into true and false datatype
                feat = textonfeat;

            end
            save(strcat(pwd,'/Basal_hist_data_at_different_window_size/cancer_vs_normal_texton_feat_around_stroma_dist_',num2str(stromawidth),'.mat'),...
                'nmtextonfeat','g3textonfeat','textonfeat','numericclass','g4textonfeat','bphtextonfeat','hgptextonfeat','-v7.3');
        end
    end
    
    nstromawidth = length(stromawidtharr);
    %--Perform cancer vs non-cancer classification--
    %Generate the cancer and non-cancer label from the class label
    kval = 10;
    nsamples = length(numericclass(:));
    for stromaidx = 1:nstromawidth
        stromawidth =  stromawidtharr(stromaidx);
        load(strcat(pwd,'/Basal_hist_data_at_different_window_size/cancer_vs_normal_texton_feat_around_stroma_dist_',num2str(stromawidth),'.mat'));
        iscancer = (numericclass==3)|(numericclass==4);
        %Search for optimal paramters for the SVM
        c = cvpartition(nsamples,'KFold',20);
        minfn = @(z)kfoldLoss(fitcsvm(textonfeat,iscancer,'CVPartition',c,...
            'KernelFunction','rbf','Standardize',true,...
            'BoxConstraint',exp(z(2)),'KernelScale',exp(z(1)))); %Use the misclassification error as the objective for parameter searching
       % options = psoptimset('UseParallel',true); % Display iterative output
        [zmin fval] = patternsearch(minfn,randn(2,1),[],[]);
        %Now, go  with the minimizer found
        svmstruct = fitcsvm(textonfeat,iscancer,'KernelFunction','rbf','Standardize',true,...
            'BoxConstraint',exp(zmin(2)),'KernelScale',exp(zmin(1)));
        yout = predict(svmstruct,textonfeat);
        %Perform cross-validatition
        CVSVMModel = crossval(svmstruct); %Perform cross-validation
        %Compute mis-classification
        misclass = kfoldLoss(CVSVMModel);
        disp(['Stroma width: ' num2str(stromawidth),', Box constrain C: ', num2str(exp(zmin(2))), ', Kernel scale: ',...
            num2str(exp(zmin(1))),', misclassification: ',num2str(misclass)])
        [cf_mat,class] = confusionmat(CVSVMModel.Y,kfoldPredict(CVSVMModel));%Compute confusion matrices on the samples that was not used for training
        cf_mat = cf_mat./repmat(sum(cf_mat,2),[1 2]);
        disp(['Confusion matrix on testing data: ']);
        cf_mat
        %Compute the confusion matrix 
        
    end
%end

