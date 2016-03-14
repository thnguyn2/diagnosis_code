%This program generate the data and perform PCA and FDL to find the
%direction with maximum variance (PCA) and that with maximum separation (FLD)
%The example comes up in 2D but you should know how to convert it to n-d
%Author: Tan H.Nguyen
%Email: thnguyn2@illinois.edu
%==================This program compute the direction for maximizing
%discrimination between 2 data groups================================
function [w,x1_proj,x2_proj]=fld(data1,data2,k)
%w: projecting direction for maximizing discrimination
%x1_proj: projection of data1 to u
%x2_proj: projection of data2 to u
%k: dimension of the low-dimensional space that we are projecting to

%Within class variance
s1=cov(data1');
s2=cov(data2');
sw=s1+s2;
%Compute between class scatter
data=[data1,data2];

%Average of all data
mu=mean(data,2);
mu1=mean(data1,2);
mu2=mean(data2,2);
n1=size(data1,2);
n2=size(data2,2);
sb=n1*(mu1-mu)*(mu1-mu)'+n2*(mu2-mu)*(mu2-mu)';
%Solve for the w in the generalize eigen value problem and just keep the
%best direction
[w,value]=eig(sb,sw);

w = w(:,2:k+1);
%Plot the direction that maximize class discrimination
%Find projection of data on to w
x1_proj=w'*data1;
x2_proj=w'*data2;
%Plot projected data
end
