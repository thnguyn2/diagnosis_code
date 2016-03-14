clc;
clear all;
close all;
cf_dir= './confusionmatrices/'

%This code compute the interpolated values for smaller sample size
xs_interp = linspace(0,1,200);
addpath(cf_dir);
load('23cf.mat');
uniq = [true;diff(xs)~=0];
ys_interp_23 = interp1(xs(uniq),ys(uniq),xs_interp,'linear');

figure(1);
plot(xs_interp,ys_interp_23,'--r');hold on;
auc_arr=zeros(16,1);
auc_arr(1)=auc;
s1 = sprintf('Gleason 2+3, AUC = %1.2f',auc);

load('32cf.mat');
uniq = [true;diff(xs)~=0];
ys_interp_32 = interp1(xs(uniq),ys(uniq),xs_interp,'linear');

plot(xs_interp,ys_interp_32,'--b');hold on;
auc_arr(2)=auc;
s2 = sprintf('Gleason 3+2, AUC = %1.2f',auc);


load('33cf.mat');
uniq = [true;diff(xs)~=0];
ys_interp_33 = interp1(xs(uniq),ys(uniq),xs_interp,'linear');

plot(xs_interp,ys_interp_33,'--k');hold on;
auc_arr(3)=auc;
s3 = sprintf('Gleason 3+3, AUC = %1.2f',auc);


load('34cf.mat');
uniq = [true;diff(xs)~=0];
ys_interp_34 = interp1(xs(uniq),ys(uniq),xs_interp,'linear');

plot(xs_interp,ys_interp_34,'--g');hold on;
auc_arr(4)=auc;
s4 = sprintf('Gleason 3+4, AUC = %1.2f',auc);

load('43cf.mat');
uniq = [true;diff(xs)~=0];
ys_interp_43 = interp1(xs(uniq),ys(uniq),xs_interp,'linear');

plot(xs_interp,ys_interp_43,'r');hold on;
auc_arr(5)=auc;
s5 = sprintf('Gleason 4+3, AUC = %1.2f',auc);

load('44cf.mat');
uniq = [true;diff(xs)~=0];
ys_interp_44 = interp1(xs(uniq),ys(uniq),xs_interp,'linear');

plot(xs_interp,ys_interp_44,'b');hold on;
auc_arr(6)=auc;
s6 = sprintf('Gleason 4+4, AUC = %1.2f',auc);

load('45cf.mat');
uniq = [true;diff(xs)~=0];
ys_interp_45 = interp1(xs(uniq),ys(uniq),xs_interp,'linear');

plot(xs_interp,ys_interp_45,'k');hold on;
auc_arr(7)=auc;
s7 = sprintf('Gleason 4+5, AUC = %1.2f',auc);

load('53cf.mat');
uniq = [true;diff(xs)~=0];
ys_interp_53 = interp1(xs(uniq),ys(uniq),xs_interp,'linear');

plot(xs_interp,ys_interp_53,'g');hold on;
auc_arr(8)=auc;
s8 = sprintf('Gleason 5+3, AUC = %1.2f',auc);

load('54cf.mat');
uniq = [true;diff(xs)~=0];
ys_interp_54 = interp1(xs(uniq),ys(uniq),xs_interp,'linear');

plot(xs_interp,ys_interp_54,'-.r');hold on;
auc_arr(9)=auc;
s9 = sprintf('Gleason 5+4, AUC = %1.2f',auc);

load('55cf.mat');
uniq = [true;diff(xs)~=0];
ys_interp_55 = interp1(xs(uniq),ys(uniq),xs_interp,'linear');

plot(xs_interp,ys_interp_55,'-.b');hold on;
auc_arr(10)=auc;
s10 = sprintf('Gleason 5+5, AUC = %1.2f',auc);

load('BPcf.mat');
uniq = [true;diff(xs)~=0];
ys_interp_bp = interp1(xs(uniq),ys(uniq),xs_interp,'linear');

plot(xs_interp,ys_interp_bp,':r');hold on;
auc_arr(11)=auc;
s11 = sprintf('BPH, AUC = %1.2f',auc);

load('NMcf.mat');
uniq = [true;diff(xs)~=0];
ys_interp_nm = interp1(xs(uniq),ys(uniq),xs_interp,'linear');

plot(xs_interp,ys_interp_nm,':b');hold on;
auc_arr(12)=auc;
s12 = sprintf('Normal, AUC = %1.2f',auc);

load('HGcf.mat');
uniq = [true;diff(xs)~=0];
ys_interp_hg = interp1(xs(uniq),ys(uniq),xs_interp,'linear');

plot(xs_interp,ys_interp_hg,':k');hold on;
auc_arr(12)=auc;
s13 = sprintf('HGPIN, AUC = %1.2f',auc);
h=xlabel('P_F(stroma)');
set(h,'FontSize',14);
h=ylabel('P_D(gland)');
set(h,'FontSize',14);
h=legend(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13);
set(h,'FontSize',14);
set(gca,'FontSize',12);%Save the data into an Excel sheet
temp_mat = [xs_interp(:) ys_interp_23(:) ys_interp_32(:) ys_interp_33(:) ys_interp_34(:) ys_interp_43(:) ys_interp_44(:) ...
    ys_interp_45(:) ys_interp_53(:) ys_interp_54(:) ys_interp_55(:) ys_interp_bp(:) ys_interp_nm(:) ys_interp_hg(:)];
xlswrite('.\confusionmatrices\rocs.xls', temp_mat)
