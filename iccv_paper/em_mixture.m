function [cpi,calpha,gmu,gstd,Ncl,L]=em_mixture(Di,Fi,L)
    %implements the EM algorithm for diagnosis
    %Author: Tan H. Nguyen
    %Last update:
    %Inputs:
    %   Y: an observation. Matrix of n x d elements
    %   Di:  a N x 1 vector where each entry is an index of the sample in
    %   the score list
    %   Fi: a cell array. Each array element is a Ni[k] x P cell elements
    %   that contains the features of the glands. P is the number of features
    %   for each gland
      
    %Outputs:
    %   gmu, gstd: P x L matrices for the mean and standard
    %   deviation of L classes for the glands
    %   calpha: a Ncl x 1 vector of the prior for Ncl different diagnosis
    %   outcmes. Ncl = length(unique([X1i, X2i],'rows')
    %   cpi: a Ncl x 1 vector of the prior for all the diagnosis outcomes
    %   L: number of classes of the invidiual glands 
    %   calp: a Ncl x1 vector for the mixture for each diagnosis outcome.
    %   Ncl: number of diagnosis outcomes in the training set
    P = size(Fi{1},2);%Number of features
    gmu = zeros(L,P); %Mean and standard deviation for each class of the gland
    gstd = zeros(L,P);
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
    gmu = randn(P,L);
    gstd = ones(P,L);
    alphak = alpha;
    gmuk = gmu;
    gstdk = gstd;
    %Display information on the estimation before running the EM algorithm
    disp('Mixture portion: ')
    disp(num2str(alphak));
    disp('Mean: ')
    disp(num2str(gmuk));
    disp('Std: ')
    disp(num2str(gstdk));
    niter = 100;
    llarray = zeros(1,0);
    for iter=1:niter
        %First, update the prior distribution for the Gleason score
        pik1 = zeros(D,1);
        for didx=1:D
            %Find all the sample that has current diagnosis outcome d
            matchsampleidx=find(Di==didx);
            ngoodsample = length(matchsampleidx);
            for matchidx=1:ngoodsample %Sum from 1 to N
                for lidx=1:L
                    tempp = compute_p1_given_fi_d(lidx,didx,Fi{matchsampleidx(matchidx)},alphak,gmuk,gstdk);
                    pik1(didx) = pik1(didx) + sum(tempp); %The second sum summs over ik
                end
            end
        end
        pik1 = pik1/sum(pik1);

        %Next, update the parameter for alpha matrix
        alphak1 = zeros(D,L);
        for didx=1:D
            matchsampleidx=find(Di==didx);
            ngoodsample = length(matchsampleidx);
            for lidx=1:L
                for matchidx=1:ngoodsample %Sum from 1 to N
                    tempp = compute_p1_given_fi_d(lidx,didx,Fi{matchsampleidx(matchidx)},alphak,gmuk,gstdk);
                    alphak1(didx,lidx)=alphak1(didx,lidx) + sum(tempp);        
                end

            end
        end 
        alphak1 = alphak1./repmat(sum(alphak1,2),[1 L]);

        %Update the parameters for each gland
        gmuk1 = zeros(P,L);

        for lidx=1:L
            num = zeros(P,1);
            den = 0;
            for didx=1:D
                matchsampleidx=find(Di==didx);
                ngoodsample = length(matchsampleidx);
                for matchidx=1:ngoodsample %Sum from 1 to N
                    tempp = compute_p1_given_fi_d(lidx,didx,Fi{matchsampleidx(matchidx)},alphak,gmuk,gstdk);
                    fp = Fi{matchsampleidx(matchidx)}.*repmat(tempp,[1 P]);
                    tempval = sum(fp,1);
                    num = num + tempval';
                    den = den + sum(tempp);
                end
            end
            gmuk1(:,lidx) = num/den;
        end

        %Next, update the parameters for the standard deviation
        gstdk1 = zeros(P,L);

        for lidx=1:L
            num = zeros(P,1);
            den = 0;
            for didx=1:D
                matchsampleidx=find(Di==didx);
                ngoodsample = length(matchsampleidx);
                for matchidx=1:ngoodsample %Sum from 1 to N
                    tempp = compute_p1_given_fi_d(lidx,didx,Fi{matchsampleidx(matchidx)},alphak,gmuk,gstdk);
                    Ni = size(Fi{matchsampleidx(matchidx)},1);
                    curmean = gmuk1(:,lidx);curmean = curmean(:)';
                    fp = ((Fi{matchsampleidx(matchidx)}-repmat(curmean,[Ni 1])).^2).*repmat(tempp,[1 P]);
                    tempval = sum(fp,1);
                    num = num + tempval';
                    den = den + sum(tempp);
                end
            end
            gstdk1(:,lidx) = sqrt(num/den);
        end
        alphak = alphak1;
        pik = pik1;
        gmuk = gmuk1;
        gstdk = gstdk1;
        disp(['Current iteration: ' num2str(iter)]);
        %disp('Outcome prior: ')
        %disp(num2str(pik));
        disp('Mixture portion: ')
        disp(num2str(alphak));
        disp('Mean: ')
        disp(num2str(gmuk));
        disp('Std: ')
        disp(num2str(gstdk));
        

    end
    
end

function pld = compute_p1_given_fi_d(l,d,Fi,alphak,gmuk,gstdk)
    %Compute p_thetak(G_ik = l|Di,fik): the likelihood of the gland ik
    %of the sample i have label l-th given its diagnosis results Di
    %and its feature fik
    %Inputs:
    %   Fi: a Ni x P: a matrix of Ni glands corresponding to a the class
    %   i-th
    %   d: a scalar index that specifies the outcome of the core that we wants to compute the likelihood 
    %   l: a scalar value indicating the class of the gland that we are
    %   computing, which specifies the column of alphak, gmuk, gstdk
    %   gmuk, gstdk:  L x P arraies of the parameters for L different
    %   classes and P features
    %   alphak: an D x L vector for the likelihood of all possible
    %   outcome of a core
    %Outputs:
    % a Ni x 1 matrix where each entry is the likelihood for sample ik (from 1 to Ni) i.e. p_thetak(G_ik = l|D,fik)
    Ni = size(Fi,1); %Number of glands in the current core
    P = size(Fi,2); %Mumber of features
    L = size(gmuk,2);
    p_arr = zeros(Ni,L); %This matrix contains p(fik|l) for l=1:L and ik = 1:Nk
    for lp = 1:L
        curmu = gmuk(:,lp);
        curstd = gstdk(:,lp);
        tempval = (Fi'-repmat(curmu,[1,Ni]));
        tempval = tempval./repmat(curstd,[1,Ni]);
        tempval = tempval.^2;
        tempval = tempval/2;
        tempval = sum(tempval,1);
        p_arr(:,lp) = (exp(-tempval'))/(2*pi)^(P/2)/det(diag(curstd));
    end
    p_alpha = p_arr.*repmat(alphak(d,:),[Ni 1]);
    pld = p_alpha(:,l)./sum(p_alpha,2);
end

%The following function is not yet tested.....Not working properly yet
function [lnpy] = compute_likelihood_yi(Di,Fi,alphak,gmuk,gstdk,pik)
%This function compute the likelihood of the data P_theta^k(di,fi)
%Inputs:
%   Di: a N x 1 vector of the diagnosis outcome for each core
%   Fi: a N x 1 cell vector where each cell is a Ni x P matrix containings
%   all the features. Here Ni is the number of glands and P is the number
%   of features of that core
%   alphak, gmuk, gstdk: see definition above
%   pik: prior p(di)
%Outputs:
%   lnpy: a scalar vector corresponding to the log likelihood of the
%   observation
%Outputs:
    N = size(Fi,1); %Number of cores
    L = size(alphak,2); %Number of sub classes
    finres = 0;
    for nidx=1:N
        curF = Fi{nidx};%Get the current feature
        di = Di(nidx);
        pdi = pik(di);
        Ni = size(curF,1);
        P = size(curF,2);
        tempres = zeros(1,L); %This matrix store p(fik,di,l)
        for lidx=1:L
           pl_di = alphak(di,lidx);
           
           %This code is for computing p_theta_k(fi|l)
           curmean = gmuk(:,lidx);
           curstd = gstdk(:,lidx);
           errorvect = curF -repmat(curmean(:)',[Ni,1]);
           errorvect=errorvect./repmat(curstd(:)',[Ni,1]);
           errorvect = (errorvect.^2)/2;
           errorvect = sum(errorvect,2);%This sum is done over the features
           pfi_l = exp(-errorvect)/(2*pi)^(P/2)/det(diag(curstd));
          
           inprod = exp(sum(log(pfi_l))); %Compute the product of p_theta_k(fik|l) for different glands
           tempres(:,lidx) = inprod*pl_di*pdi;
        end
        %Summing over the class L and take the log
        tempres = sum(tempres);
        
        %Summing over the samples
        finres = finres + log(tempres);
    end
    lnpy = finres;
end

