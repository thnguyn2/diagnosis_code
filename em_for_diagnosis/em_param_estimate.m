function [cpi,calpha,gmu,gcov]=em_param_estimate(Di,Fi,L,niter)
    %implements the EM algorithm for diagnosis
    %Author: Tan H. Nguyen
    %Last update:
    %Inputs:
    %   Di:  a N x 1 vector where each entry is an index of the sample in
    %   the score list
    %   Fi: a cell array. Each array element is a Ni[k] x P cell elements
    %   that contains the features of the glands. P is the number of features
    %   for each gland
    %   niter: number of iterations for the EM
    %   L: number of classes of the invidiual glands 
    %Outputs:
    %   gmu, gstd: P x L matrices for the mean and standard
    %   deviation of L classes for the glands
    %   calpha: a Ncl x 1 vector of the prior for Ncl different diagnosis
    %   outcmes. Ncl = length(unique([X1i, X2i],'rows')
    %   cpi: a Ncl x 1 vector of the prior for all the diagnosis outcomes

    %   calp: a Ncl x1 vector for the mixture for each diagnosis outcome.
  
    P = size(Fi{1},2);%Number of features
    gmu = cell(L,1); %Mean and standard deviation for each class of the gland
    gcov = cell(L,1);
   
    N = size(Fi,1);
    D = length(unique(Di));
    P = size(Fi{1},2);
    %First, get the number of glands per core
    Ni = zeros(N,1); %Number of glands per core
    for coreidx=1:N
        Ni(coreidx)=size(Fi{coreidx},1);        
    end
 
    %Initialize the parameter
    alpha = 0.5*ones(D,L);
    
    for lidx=1:L
        gmu{lidx}=randn(P,1);
        gcov{lidx}=eye(P,P); %Initialize the covariance matrix with the identity matrix
    end
    alphak = alpha;
    gmuk = gmu;
    gcovk = gcov;
    
    for iter=1:niter
        %First, update the prior distribution for the Gleason score
        pik1 = zeros(D,1);
        for didx=1:D
            %Find all the sample that has current diagnosis outcome d
            matchsampleidx=find(Di==didx);
            ngoodsample = length(matchsampleidx);
            pik1(didx)=ngoodsample;
        end
        pik1 = pik1/sum(pik1);
        %Next, update the parameter for alpha matrix
        alphak1 = zeros(D,L);
        for didx=1:D
            matchsampleidx=find(Di==didx);
            ngoodsample = length(matchsampleidx);
            for lidx=1:L
                for matchidx=1:ngoodsample %Sum from 1 to N
                    tempp = compute_p1_given_fi_d(lidx,didx,Fi{matchsampleidx(matchidx)},alphak,gmuk,gcovk);
                    alphak1(didx,lidx)=alphak1(didx,lidx) + sum(tempp);        
                end

            end
        end 
        alphak1 = alphak1./repmat(sum(alphak1,2),[1 L]);

        %Update the parameters for each gland
        gmuk1 = cell(L,1);

        for lidx=1:L
            num = zeros(P,1);
            den = 0;
            for didx=1:D
                matchsampleidx=find(Di==didx);
                ngoodsample = length(matchsampleidx);
                for matchidx=1:ngoodsample %Sum from 1 to N
                    tempp = compute_p1_given_fi_d(lidx,didx,Fi{matchsampleidx(matchidx)},alphak,gmuk,gcovk);
                    fp = Fi{matchsampleidx(matchidx)}.*repmat(tempp,[1 P]);
                    tempval = sum(fp,1);
                    num = num + tempval';
                    den = den + sum(tempp);
                end
            end
            gmuk1{lidx} = num/den;
        end

        %Next, update the parameters for the standard deviation
        gcovk1 = cell(L,1);

        for lidx=1:L
            num = zeros(P,P);
            den = 0;
            for didx=1:D
                matchsampleidx=find(Di==didx);
                ngoodsample = length(matchsampleidx);
                for matchidx=1:ngoodsample %Sum from 1 to N
                    plgfi = compute_p1_given_fi_d(lidx,didx,Fi{matchsampleidx(matchidx)},alphak,gmuk,gcovk);
                    Ni = size(Fi{matchsampleidx(matchidx)},1);
                    curmean = gmuk1{lidx};
                    curmean = curmean(:)';
                    tempdiff = Fi{matchsampleidx(matchidx)}-repmat(curmean,[Ni 1]);
                    for glandidx=1:Ni
                       curerr = tempdiff(glandidx,:);
                       curerr = curerr(:); %fik-mu_k1
                       num = num + plgfi(glandidx)*(curerr*curerr');
                    end
                    den = den + sum(plgfi);
                end
            end
            gcovk1{lidx} = num/den;
        end
        alphak = alphak1;
        pik = pik1;
        gmuk = gmuk1;
        gcovk = gcovk1;
        disp(['Current iteration: ' num2str(iter)]);


    end

    %Return the set of parameter
    cpi = pik;
    calpha = alphak;
    gmu = gmuk;
    gcov = gcovk;
    disp('Outcome prior: ')
    disp(num2str(pik));
    disp('Mixture portion: ')
    disp(num2str(alphak));
    disp('Mean 1')
    disp(num2str(gmuk{1}));
    disp('Mean 2')
    disp(num2str(gmuk{2}));
    disp('cov  1')
    disp(num2str(gcovk{1}));
    disp('cov  2')
    disp(num2str(gcovk{2}));
        
end

