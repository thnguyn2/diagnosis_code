function [pstroma,pgland,nstroma,ngland,binarr]=computelikelihood(fim,lblim,...
    minval,maxval,nbins)
%Compute the likelihood p(fim|x=stroma) and p(fim|x=gland)
%Inputs:
%   fim: functional margin image
%   lblim: manually labelled image
%   maxval, minval: maximum and mininum values of the histogram bins
%   nbins: number of histogram bins
%Outputs:
%   pstroma, pglands: histogram distribution of glands and stromas
%   nstroma, nglands: number of stroma and gland pixels
%------------------------------------------------------------------
    stromaidx = find(lblim==-1);
    glandidx = find(lblim==1);
    nstroma = numel(stromaidx);
    ngland = numel(glandidx);
    binarr = linspace(minval,maxval,nbins);
    pstroma = hist(fim(stromaidx),binarr)/nstroma;
    pgland = hist(fim(glandidx),binarr)/ngland;
end