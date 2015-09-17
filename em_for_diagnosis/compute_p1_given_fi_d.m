function pld = compute_p1_given_fi_d(l,d,Fi,alphak,gmuk,gcovk)
    %Compute p_thetak(G_ik = l|Di,fik): the likelihood of the gland ik
    %of the sample i have label l-th given its diagnosis results Di
    %and its feature fik
    %Inputs:
    %   Fi: a Ni x P: a matrix of Ni glands corresponding to a the class
    %   i-th
    %   d: a scalar index that specifies the outcome of the core that we wants to compute the likelihood 
    %   l: a scalar value indicating the class of the gland that we are
    %   computing, which specifies the column of alphak, gmuk, gstdk
    %   gmuk:  L x 1 vector of cells of the parameters for L different
    %   classes. Each cell is a P x 1 matrix of the mean. This is the mean
    %   for all the grades.
    %   gcovk: an L x 1 vector of cells where each cell is a correlation
    %   matrix of size P x P
    %   classes and P features
    %   alphak: an D x L vector for the likelihood of all possible
    %   outcome of a core
    %Outputs:
    % a Ni x 1 matrix where each entry is the likelihood for sample ik (from 1 to Ni) i.e. p_thetak(G_ik = l|D,fik)
    Ni = size(Fi,1); %Number of glands in the current core
    P = size(Fi,2); %Mumber of features
    L = size(gmuk,1);%Number of classes
    p_arr = zeros(Ni,L); %This matrix contains p(fik|l) for l=1:L and ik = 1:Nk
    expo_term = zeros(Ni,1);%The expoential term for each glands
    scalingfact= ((2*pi)^(-P/2));
    for lp = 1:L
        curmu = gmuk{lp};
        curcov = gcovk{lp};
        invcurcov = inv(curcov);
        detsqrt2=(det(curcov)^(-0.5));
        tempval = (Fi'-repmat(curmu,[1,Ni]));      
        for ik=1:Ni
            curvect = tempval(:,ik);
            expo_term(ik)=-0.5*curvect'*invcurcov*curvect;
        end
        p_arr(:,lp) = scalingfact*detsqrt2*exp(expo_term);
    end
    p_alpha = p_arr.*repmat(alphak(d,:),[Ni 1]);
    pld = p_alpha(:,l)./sum(p_alpha,2);
end
