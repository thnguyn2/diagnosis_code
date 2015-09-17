%This program implements the EM algorithm for the cancer diagnosis
function demo
    clc;
    close all;
    clear all;
    Ncl = 100; %Number of sample per diagnosis results
    alpha = [1 0;%Each row is for a diagnosis result, each column is for a gland's grade
            0 1
            ];
   mu = cell(size(alpha,2),1);
   mu{1}=[2;1]; mu{2}=[1.8;1.2];
   cov = cell(size(alpha,2),1);
   cov{1}=[8  0.2; %Each column is a feature
            0.2 5];
   cov{2}=[8 0.5;
           0.5 4];
   
   [Di,Ni,Fi,N]=generate_observation_data(Ncl,alpha,mu,cov);
   niter = 30;
   %The following two parameters should be chosen ahead
   L = 2; %Number of possible grades for the glands
   D = length(unique(Di)); %Number of possible diagnosis outcomes
   [cpi,calpha,gmu,gcov]= em_param_estimate(Di,Fi,L,niter);
   c = 1; %Diagnosis result of interest
   N = size(Fi,1);%Number of samples
   pc_g_f = zeros(N,D);%p(c|f) for all N samples and D possible outcomes
   for sampleidx=1:N
    pc_g_f(sampleidx,:) = compute_diagnosis_likelihood(Fi{sampleidx},cpi,calpha,gmu,gcov);
   end
   [maxval,diag_res] = max(pc_g_f,[],2);
   confmat = confusionmat(Di,diag_res);
   disp('Confusion mat: ');
   disp(num2str(confmat));
end

function pc_g_f = compute_diagnosis_likelihood(Fi,cpi,calpha,gmu,gcov)
    %Compute the likelihood of the diagnosis results given the core
    %Inputs:
    %   Fi: a Ni[k] x P cell element of the feature. Each Fi is the feature
    %   for a core.
    %   that contains the features of the glands. P is the number of features
    %   for each gland
    %   gmu, gstd: P x L matrices for the mean and standard
    %   deviation of L classes for the glands
    %   calpha: a Ncl x 1 vector of the prior for Ncl different diagnosis
    %   outcmes. Ncl = length(unique([X1i, X2i],'rows')
    %   cpi: a Ncl x 1 vector of the prior for all the diagnosis outcomes
    %   d: the index for the diagnosis outcome that we area considering
    %Outputs:
    %   pc_g_f: a D x 1 matrix for the p(c|Fi), the likelihood of the
    %   diagnosis results given the feature Fi
    L = size(calpha,2); %Number of possible grades for the glands
    D = size(calpha,1); %Number of possible diagnosis results
    P = size(Fi,2); %Number of features
    Ni = size(Fi,1); %Number of glands for each core
    scalingfact= ((2*pi)^(-P/2));
   
    p_fi_g_c = zeros(D,1); %p(Fik|c) likelihood of the feature given the diagnosis results c for different values of c
    for cidx=1:D
        p_fik_l_g_ci = zeros(Ni,L);%This is p(fik,gi|cidx)
        for lidx=1:L
            p_l_g_c =  calpha(cidx,lidx);
            p_fik_g_l = zeros(Ni,1);%P(Fik|l). Each row is for 1 gland
            %Compute p(Fi|l) for all different values of L
             curmu = gmu{lidx};
             curcov = gcov{lidx};
             invcurcov = inv(curcov);
             detsqrt2=(det(curcov)^(-0.5));
             for glandidx=1:Ni
                curfeat = Fi(glandidx,:);
                curfeat = curfeat(:);
                curvect = curfeat - curmu;
                expoterm = -0.5*curvect'*invcurcov*curvect;
                p_fik_g_l(glandidx,1) = scalingfact.*detsqrt2*exp(expoterm);
             end
            p_fik_l_g_ci(:,lidx) =   p_l_g_c*p_fik_g_l;
        end
        p_fik_g_c = sum(p_fik_l_g_ci,2);%p(Fik|c) - summing over all possible class of a core.
        p_fi_g_c(cidx) = exp(sum(log(p_fik_g_c)));
    end
    p_fi_c =  p_fi_g_c(:).*cpi(:);%p(Fi,c)
    den = sum(p_fi_c);
    pc_g_f=p_fi_c/den;
    
end


