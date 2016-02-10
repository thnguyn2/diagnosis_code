function [cpi,calpha,gmu,gstd,Ncl,L]=em_mixture(Di,Fi,L)
    %implements the EM algorithm for diagnosis. The data is assume to be a
    %mixture of exponential distribution
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
    if (P>=2)
        error('Not supported');
    end
    lambda = rand(1,L); %Mean and standard deviation for each class of the gland
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
    alphak = alpha;
    lambdak = lambda;
    %Display information on the estimation before running the EM algorithm
    disp('Mixture portion: ')
    disp(num2str(alphak));
    disp('Lambda: ')
    disp(num2str(lambdak));
    
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
                    tempp = compute_p1_given_fi_d(lidx,didx,Fi{matchsampleidx(matchidx)},alphak,lambdak);
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
                    tempp = compute_p1_given_fi_d(lidx,didx,Fi{matchsampleidx(matchidx)},alphak,lambdak);
                    alphak1(didx,lidx)=alphak1(didx,lidx) + sum(tempp);        
                end

            end
        end 
        alphak1 = alphak1./repmat(sum(alphak1,2),[1 L]);

        %Update the parameters for each gland
        lambdak1 = zeros(1,L);

        for lidx=1:L
            num = zeros(P,1);
            den = 0;
            for didx=1:D%D is the number of diagnosis outcomes
                matchsampleidx=find(Di==didx);
                ngoodsample = length(matchsampleidx);
                for matchidx=1:ngoodsample %Sum from 1 to N
                    tempp = compute_p1_given_fi_d(lidx,didx,Fi{matchsampleidx(matchidx)},alphak,lambdak);
                    fp = Fi{matchsampleidx(matchidx)}.*repmat(tempp,[1 P]);
                    tempval = sum(fp,1);%This sum is over the glands in each core
                    num = num + sum(tempp);
                    den = den + tempval;
                end
            end
            lambdak1(:,lidx) = num/den;
        end

        alphak = alphak1;
        pik = pik1;
        lambdak = lambdak1;
        disp(['Current iteration: ' num2str(iter)]);
        %disp('Outcome prior: ')
        %disp(num2str(pik));
        disp('Mixture portion: ')
        disp(num2str(alphak));
        disp('Mean of the dist: ')
        disp(num2str(1./lambdak));
        

    end
    
end

function pld = compute_p1_given_fi_d(l,d,Fi,alphak,lambdak)
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
    L = length(lambdak);
    p_arr = zeros(Ni,L); %This matrix contains p(fik|l) for l=1:L and ik = 1:Nk
    for lp = 1:L %Go through different possible glands values
        curlambda = lambdak(:,lp);
%         tempval = (Fi'-repmat(curmu,[1,Ni]));
%         tempval = tempval./repmat(curstd,[1,Ni]);
%         tempval = tempval.^2;
%         tempval = tempval/2;
%         tempval = sum(tempval,1);
%         p_arr(:,lp) = (exp(-tempval'))/(2*pi)^(P/2)/det(diag(curstd));
        p_arr(:,lp)=curlambda*exp(-curlambda*Fi(:));
    end
    p_alpha = p_arr.*repmat(alphak(d,:),[Ni 1]);
    pld = p_alpha(:,l)./sum(p_alpha,2);
end
