%This program implements the EM algorithm for the cancer diagnosis
function demo
    clc;
    close all;
    clear all;
   %[Di,Ni,Fi,N]=generate_observation_data();
   %First, read the diagnosis data
   
    addpath(strcat(cd(cd('..')),'\support'));
    datapath = 'F:\TMA_cores_and_diagnosis\';
    diagpath = 'F:\TMA_cores_and_diagnosis\diag\';
    featpath = 'F:\TMA_cores_and_diagnosis\feat\';%This is the folder for analysing diagnosis features
    disp('Finding the list of all labels');
    %[filenames,glandnames,gradelist]=findFileNameFromROIs(datapath);
     load glandmorpfeat.mat;
     %Show the histogram of features based on the Gleason grades     
     [Di,Fi]=compute_feature_statistics_cancer_class_for_em(glandmorpfeat,classvect,true,indivglandfeat);%Draw the histogram of features for the cancer group
    %Draw the class distribution
  
   % draw_data(Di,Fi);
     %Find the grade 
    em_mixture(Di,Fi,2); %Mixtures of Gaussian
    
    %Mixture of expoential
    %em_mixture_exp_dist(Di,Fi,2)
    
end

%Draw the datra of different classes
function draw_data(Di,Fi)
    nbins = 50;
    for featidx=1
        titleval = strcat(['Feat ' num2str(featidx)]);
        title(titleval);
        for classidx=1:4
            figure(featidx);
       
            subplot(2,2,classidx);
            sampleidx = find(Di==classidx);
            ni = length(sampleidx);
            minval = 1e+6;
            maxval = -1e+6;
            %Find the min and max of current features to compute the
            %histogram
            for tempidx=1:ni
                curfeat = Fi{sampleidx(tempidx)}(:,featidx);
                curmin = mean(curfeat);
                curmax = max(curfeat);
                if (curmin<minval)
                    minval = curmin;
                end
                if (curmax>maxval)
                    maxval = curmax;
                end               
            end
            %Generate the histogram bins
            binval = linspace(minval,maxval*0.7,nbins);
            for tempidx=1:ni
                curfeat = Fi{sampleidx(tempidx)}(:,featidx);
                [curhist] = hist(curfeat,binval);
                curhist = curhist/sum(curhist);   
                plot(binval,curhist);
                hold on;
                drawnow;
            end
            
        end
    end
    

end
