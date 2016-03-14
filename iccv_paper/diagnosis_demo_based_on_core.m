%function diagnosis_demo_based_on_core
%This function perform the analysis of features to discriminate Gleason
%grade 3 and 4 with the groundtruth generated by Vicky. The features to run
%are based on the labeled image. We will mostly evaluate the fusion and
%cribriform pattern only. 
%Author: Tan H. Nguyen
%Last update: March 8th, 2015
%University of Illinois at Urbana-Chamapaign
    clc;
    clear all;
    close all;
    datafolder = 'E:\TMA_cores_and_diagnosis\diagnosis_of_vicky\';
    addpath('..\support');
    g3folder =strcat(datafolder,'g3\');
    g4folder = strcat(datafolder,'g4\');
    g3files = dir(strcat(g3folder,'*_lbl_*g3.tif'));
    g4files = dir(strcat(g4folder,'*_lbl_*g4.tif'));
    n3 = size(g3files,1); %Number of G3 samples
    n4 = size(g4files,1); %Number of G4 samples
    param.smallestlumenarea=2000;
    param.minglandarea = 5000;
    
    param.glswstd = 25; %Ls and g windows size
    param.stromawidth = 15;
    param.cleftingwidth =90;
    nfeattype = 8;
    param.nbins = 3;
    param.nfeatwithhist=2;
    ntextons = 50;
    g3feat = zeros(0,nfeattype);
    g4feat = zeros(0,nfeattype);
    g3textonfeat=zeros(0,ntextons);
    g4textonfeat=zeros(0,ntextons);
    iter = 0 ;
     tic;
    for idx=1:max(n3,n4)
        iter = iter+1;
        %Step 1: count the number of true glands and the number of lumen
        %each gland has
        %First, load the label map and cut out small fused regions. Make
        %sure that a slight fusion will not be detected as fused glands.
        %This seems to have been done in the code that generate the label
        %map.
        if (idx<=n3)
            %g3files(idx).name = 'B9_lbl_1_g3.tif';
            curg3filename = strcat(g3folder,g3files(idx).name);
            tic;
            res=process_single_core(curg3filename,param);
            averageTime = toc/1000;
            disp(['Average extracting time' num2str(averageTime)]);
            disp(['G3 - ' g3files(idx).name ': ' num2str(res.mean_nlum) ', ' num2str(res.mean_dist)]);
            slashpos = strfind(g3files(idx).name,'_');
            if ((isnan(res.mean_dist))|(res.nglands<3)|(res.roisize<0.5*1e+6))
                 movefile(curg3filename,strcat(g3folder,'confused core\',g3files(idx).name));%Move the file to the new new location
                
            else
                res.mean_dist;y=res.mean_nlum;
                if (res.mean_dist>4)
                    res.mean_dist=4;
                end
                figure(2);plot(res.ratio_area,res.mean_dist,'^b','linewidth',3);hold on;
                text(res.ratio_area+0.0001,res.mean_dist,strcat(g3files(idx).name(1:slashpos(1)-1),'-',...
                    g3files(idx).name(slashpos(2)+1:slashpos(3)-1),sprintf('-%0.1f',res.roisize/1e+6)),'Color','b');
                g3feat(end+1,:)=[res.mean_dist,res.mean_cir,res.mean_nlum,res.max_nlum,res.fused_ratio,res.med_dist,res.ratio_area,0];
                %g3textonfeat(end+1,:)=res.hist_tex_idx(:)';
            end
        end
        
        if (idx<=n4)
            %g4files(idx).name = 'M4_lbl_2_g4.tif';
            curg4filename = strcat(g4folder,g4files(idx).name);
            
            res=process_single_core(curg4filename,param);
             disp(['G4 - ' g4files(idx).name ': ' num2str(res.mean_nlum) ', ' num2str(res.mean_dist)]);
            if ((isnan(res.mean_dist))|(res.nglands<3)|(res.roisize<0.5*1e+6))
                        %error(strcat(curnmfilename,'has NaN histogram'));
                        movefile(curg4filename,strcat(g4folder,'confused core\',g4files(idx).name));%Move the file to the new new location
            else
                slashpos = strfind(g4files(idx).name,'_');
                if (res.mean_dist>4)
                    res.mean_dist=4;
                end
                if (res.nglands>=1) %Do not goes with 2 here since the whole core can be a strongly fused gland
                    figure(2);plot(res.ratio_area,res.mean_dist,'or','linewidth',3);hold on;
                    text(res.ratio_area+0.0001,res.mean_dist,strcat(g4files(idx).name(1:slashpos(1)-1),'-',...
                        g4files(idx).name(slashpos(2)+1:slashpos(3)-1),sprintf('-%0.1f',res.roisize/1e+6)),'Color','r');
                    xlabel('stdarea of gland/mean area of gland');
                    ylabel('Mean distortion');
                end
                g4feat(end+1,:)=[res.mean_dist,res.mean_cir,res.mean_nlum,res.max_nlum,res.fused_ratio,res.med_dist,res.ratio_area,0];
                %g4textonfeat(end+1,:)=res.hist_tex_idx(:)';
            end
        end
        
        
            
        figure(2)
        feat = [g3feat;g4feat];
        textonfeat = [g3textonfeat;g4textonfeat];
%        feat = [feat, textonfeat];
        curn3 = size(g3feat,1);
        curn4 = size(g4feat,1);
        class = (1:(curn3+curn4))'>curn3;
        
            featcolor='mgbkcyr';
            figure(3)
            auc_arr = zeros(nfeattype,1);
            x = cell(nfeattype,1);
            y = cell(nfeattype,1);
           if (iter>20)
                for featidx=1:nfeattype
                        curfeat = feat(:,featidx);
                        b = glmfit(curfeat,class,'normal');%Generalized linear model
                        p = glmval(b,curfeat,'logit');%Compute the fitted probability
                        [x{featidx},y{featidx},t,auc_arr(featidx)] = perfcurve(class,p,true);%Compute the ROC curve by providing the groundtruth
                        plot(x{featidx},y{featidx},featcolor(mod(featidx-1,7)+1));
                        xlabel('False positive rate'); ylabel('True positive rate');
                        hold on;                
                        plot(x{featidx},y{featidx},strcat('--',featcolor(mod(featidx-1,7)+1)));


                end
                b = glmfit(feat(:,1:nfeattype),class,'normal');%Generalized linear model for all features
                p = glmval(b,feat(:,1:nfeattype),'logit');%Compute the fitted probability
                [x_c,y_c,t,auc] = perfcurve(class,p,true);%Compute the ROC curve by providing the groundtruth
                plot(x_c,y_c,'-xr','linewidth',2);
                legend(strcat('Mean dist, AUC=',num2str(auc_arr(1))),strcat('Mean circ, AUC=',num2str(auc_arr(2))),strcat('Avg. #lumen, AUC=',num2str(auc_arr(3))),...
                            strcat('Max. nlumen AUC=',num2str(auc_arr(4))),strcat('Fused ratio. nlumen AUC=',num2str(auc_arr(5))),...
                            strcat('Med. dist, AUC=',num2str(auc_arr(6))),strcat('Std area/Mean area, AUC=',num2str(auc_arr(7))),...
                            strcat('50/50 line, AUC=',num2str(auc_arr(8))),...
                            strcat('Combined, AUC=',num2str(auc)));
                hold off;
          end
    end
    averageTime = toc/1000;
    disp(['Average extracting time' num2str(averageTime)]);


        
%         %Next, try to spot out the texton that gives highest AUC for class
%         %3 and 4 and determine where to look for them
%         texton_auc_arr = zeros(ntextons,1);
%         for textonidx=1:ntextons
%                 curfeat = feat(:,end+textonidx-ntextons);
%                 b1 = glmfit(curfeat,class,'normal');%Generalized linear model
%                 p1 = glmval(b1,curfeat,'logit');%Compute the fitted probability
%                 [x1,y1,t,texton_auc_arr(textonidx)] = perfcurve(class,p1,true);%Compute the ROC curve by providing the groundtruth
%         end
        
%         figure(7);
%         subplot(121);bar(texton_auc_arr);xlabel('Texton indices');ylabel('AUC values');
%         
%        %Compute the mean and the std of two classes
%        subplot(122);
%        mean_g3 = mean(g3textonfeat,1);
%        std_g3 = std(g3textonfeat,0,1);
%        mean_g4 = mean(g4textonfeat,1);
%        std_g4 = std(g4textonfeat,0,1);
%        errorbar([1:ntextons],mean_g3,std_g3,'or');hold on;
%        errorbar([1:ntextons],mean_g4,std_g4,'^b');hold off;drawnow;
%        legend('g3','g4');
%        
%        figure(6);
%        subplot(121);
%        imagesc(g3textonfeat,[0 0.3]);title('G3 carcinoma');colorbar;title('G3 carcinoma');drawnow;
%        subplot(122);
%        imagesc(g4textonfeat,[0 0.3]);title('G4 carcinoma');colorbar;title('G4 carcinoma');drawnow;
%       
%        
%         %find 5 textons that has highest auc and display them
%         [maxval,tex_idx_sorted]=sort(texton_auc_arr,'descend');
%         texton_freq = sum(textonfeat,1)/size(textonfeat,1);
%         %subplot(122);bar(texton_freq);xlabel('Texton indices');ylabel('Frequency');
        
        
%         %Ccombine all the feature
%         b = glmfit(feat(:,1:nfeattype),class,'normal');%Generalized linear model for all features
%         p = glmval(b,feat(:,1:nfeattype),'logit');%Compute the fitted probability
%         [x_c,y_c,t,auc] = perfcurve(class,p,true);%Compute the ROC curve by providing the groundtruth
%         auc
%         plot(x_c,y_c,'-xr','linewidth',2);
%         legend(strcat('Mean dist, AUC=',num2str(auc_arr(1))),strcat('Avg. #lumen, AUC=',num2str(auc_arr(2))),...
%                     strcat('Max. nlumen AUC=',num2str(auc_arr(3))),strcat('Med. dist, AUC=',num2str(auc_arr(4))),...
%                     strcat('Mean ls, AUC=',num2str(auc_arr(5))),strcat('Mean g, AUC=',num2str(auc_arr(6))),...
%                     strcat('Combined, AUC=',num2str(auc)));
%         hold off;

%            important_loc_mask = zeros(size(res.texidx));
%            nimportfeat = 50;
%            for featidx=1:nimportfeat
%                important_loc_mask = important_loc_mask + (res.texidx==tex_idx_sorted(featidx)).*maxval(featidx);
%            end
%           
%         figure(8);
%         imagesc(res.stromastrand.*important_loc_mask); drawnow;colorbar
  
    save('roidiag.mat','g3feat','g4feat','-v7.3');
    
%end

