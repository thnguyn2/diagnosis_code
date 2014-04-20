function updateAffinityMatrix(filenames,curpath)
%This function computes the affinity matrix for each iage. This function is
%for class-separation using the FLD. With SVM classifier, we don't need
%thsi function
%Inputs:
%   filenames the name of all cores
%   curpath: path to the directory that saves response information information
%   datapath: path to the directory that saves the phase information.

%Outputs:
%   [None]

    addpath(curpath);
    
    ncols = 2048;
    nrows = 2048;
    ntextons = 80;

    h = waitbar(0,'Overal percentage...');
    k = waitbar(0,'Hough voting...');
    
    %% Computing textons for each images..... 
    nfilestotal = 0;
    for classidx=1:4 %Go through different classes
       nsamples = (length(filenames{classidx,1})); %Get the number of samples in each class
       for sampleIdx=1:nsamples
            nfilestotal =nfilestotal+1; 
            waitbar(sampleIdx/nsamples,h,'Progress...')
            cur_file_name = filenames{classidx,1}{sampleIdx,1};
            %Check to see if the label file exist
            dot_pos = strfind(cur_file_name,'.'); %Get the position of the dot
            slash_pos = strfind(cur_file_name,'\');
            label_name = cur_file_name(slash_pos(end)+1:dot_pos(1)-1);
            texton_hist_map = strcat(label_name,'_texton_hist.mat');%Name of the filter response
            affin_file = strcat(label_name,'_affin.mat');
            imname = strcat(label_name,'.mat');
            disp(['Calculating nuclei for: ' label_name ' ...']);

            if (~exist(strcat(curpath,affin_file),'file'))
                load(strcat(curpath,texton_hist_map),'histim','lblim');
                glandidx = find(lblim==1);
                stromaidx = find(lblim==2);
                lumenidx = find(lblim==0);
                ngland = length(glandidx);
                nstroma = length(stromaidx);
                nlumen = length(lumenidx);
                gland_hist = zeros(ngland,ntextons);
                stroma_hist = zeros(nstroma,ntextons);
                lumen_hist = zeros(nlumen,ntextons);
                for bandidx=1:ntextons
                    waitbar(bandidx/ntextons,k,'Progress...')
                    curband = histim(:,:,bandidx);
                    glandvect=curband(glandidx);
                    stromavect=curband(stromaidx);
                    lumenvect=curband(lumenidx);
                    gland_hist(:,bandidx)=glandvect(:);
                    stroma_hist(:,bandidx)=stromavect(:);
                    lumen_hist(:,bandidx)=lumenvect(:);
                end
                clear histim;
                
                gland_hist=gland_hist';
                lumen_hist=lumen_hist';
                stroma_hist=stroma_hist';
                mugland = sum(gland_hist,2)/ngland;
                mustroma = sum(stroma_hist,2)/nstroma;
                mulumen = sum(lumen_hist,2)/nlumen;
                mu = (mugland*ngland+mustroma*nstroma)/(ngland+nstroma);
                Sglands = (gland_hist - repmat(mugland,[1 ngland]))*(gland_hist - repmat(mugland,[1 ngland]))';
                Sstromas = (stroma_hist - repmat(mustroma,[1 nstroma]))*(stroma_hist - repmat(mustroma,[1 nstroma]))';
                Sw=Sglands+Sstromas;
                Sb = ngland*(mugland-mu)*(mugland-mu)'+nstroma*(mustroma-mu)*(mustroma-mu)';
                [V,D]=eigs(Sb,Sw,2);
                nfeatures = 1;
                W=real(V(:,1:nfeatures));
                glandproj = W'*gland_hist;
                stromaproj = W'*stroma_hist;
                if (nfeatures==1)
                    figure(4);
                    hold on;
                    plot(stromaproj','r');
                     plot(glandproj');

                end
                save(strcat(curpath,affin_file),'Sglands','Sstromas','Sw','Sb','W','gland_hist','lumen_hist','stroma_hist','-v7.3');
                %Draw the projection map between stroma & glands
                clear gland_hist;
                clear stroma_hist;
                clear lumen_hist;
                figure(1);
                imagesc(lblim);
                title('Label image');
                %Perform validation a step by showing the values
                load(strcat(curpath,texton_hist_map),'histim');
                for textonidx=1:ntextons
                    waitbar(textonidx/ntextons,k,'Progress...')
                    histim(:,:,textonidx) = histim(:,:,textonidx)*W(textonidx);
                end
                proj = sum(histim,3);
                %Mask out the lumen
                idx = find(lblim==0);
                proj(idx)=0;
                figure;
                imagesc(proj);
                colorbar;
                title('Feature value');
            
            else
                load(strcat(curpath,affin_file),'W');
                load(strcat(curpath,texton_hist_map),'histim','lblim');
                 for textonidx=1:ntextons
                    waitbar(textonidx/ntextons,k,'Progress...')
                    histim(:,:,textonidx) = histim(:,:,textonidx)*W(textonidx);
                end
                proj = sum(histim,3);
                %Mask out the lumen
                idx = find(lblim==0);
                proj(idx)=0;
                figure;
                imagesc(proj);
                colorbar;
                title(label_name);
                save(strcat(curpath,affin_file),'proj','-append');
              
            end
         
           
       end
    end
    close(k);
   
end