function [Di,Ni,Fi,N]=generate_observation_data(Ncl,alpha,mu,cov)
%Generate a dataset pf D*Ncl different cores. The dataset has D group. Each group has Ncl samples
%Inputs:
%   Ncl: number of sample for each diagnosis outcomes
%   L, D: number of classes for each core and number of diagnosis outcomes
%   alpha: a D x L likelihood matrix p(l|d). Element is the probability
%   that a gland belongs to class d comes from class l. Each row is for one
%   diagnois results, each column is for one grades of the glands
%   element is the likelihood for a gland to belongs to the major class
%   muk: a L x 1 cells where each cell is a P x 1 matrix for the mean of a
%   gland's gade
%   covk: a L x 1 cells where each cell is a P x P matrix for the cov of a
%   gland's gade

%Outputs:
%   N: the total number of samples return which is D*Ncl
%   Di:  an N x 1 array of the labels for
%   Ni: a 4Ncl x 1 array for each row is the number of glands for each
%   cores
%   Fi: a 4Ncl x 1 array. Each array element is a Ni[k] x P cell elements
%   that contains the features of the glands. P is the number of features
%   for each gland
   
    D = size(alpha,1);
    N = D*Ncl;
    
    Ni = zeros(N,1);
    Fi = cell(N,1);
    Di = zeros(N,1);
    colorarr='rgbmk';
    mu1 = mu{1}; %
    %Generate the correlation matrix
    sigma1 = cov{1};
    mu2 = mu{2};
    sigma2 = cov{2};
          
    %Compute the SVD of sigma 1 and sigma 2 so that we can generate the
    %data
    [u,s,v]=svd(sigma1);
    a1inv = u*s^0.5*u';
    [u,s,v]=svd(sigma2);
    a2inv = u*s^0.5*u';
    P = length(mu1); %Number of features
   
    for classidx=1:D
        for sampleidx=1:Ncl
            fullsampleidx=(classidx-1)*Ncl+sampleidx;
            Di(fullsampleidx)=classidx;
            Ni(fullsampleidx)=ceil(rand()*100); %Generate random number of sample
            curNi = Ni(fullsampleidx);
            curfeat = zeros(curNi,P);
            for glandidx=1:curNi
                curalpha = alpha(classidx);
                temp = rand();
                featvect = randn(P,1);
                if (temp<curalpha) %Generate the majority class
                        featvect = a1inv*featvect;
                        featvect = featvect + mu1;                    
                else %Generate sample from the minor class
                        featvect = a2inv*featvect;
                        featvect = featvect + mu2;                    
                end
                curfeat(glandidx,:)=featvect;
            end
            meanvect = mean(curfeat,1);
            if (P==2)
                figure(1);
                plot(meanvect(1),meanvect(2),strcat('+',colorarr(classidx)));hold on;
            elseif (P==1)
                plot(meanvect(1),strcat('+',colorarr(classidx)));hold on;
            
            end
            Fi{fullsampleidx}=curfeat;            
        end
    end   
end