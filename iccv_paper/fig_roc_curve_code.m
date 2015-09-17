clc;
clear all;
close all;
cf_dir= '.\confusionmatrices\'
addpath(cf_dir);
load('23cf.mat');
figure(1);
plot(xs,ys,'--r');hold on;
auc_arr=zeros(16,1);
auc_arr(1)=auc;
s1 = sprintf('Gleason 2+3, AUC = %1.2f',auc);

load('32cf.mat');
plot(xs,ys,'--b');hold on;
auc_arr(2)=auc;
s2 = sprintf('Gleason 3+2, AUC = %1.2f',auc);

load('33cf.mat');
plot(xs,ys,'--k');hold on;
auc_arr(3)=auc;
s3 = sprintf('Gleason 3+3, AUC = %1.2f',auc);

load('34cf.mat');
plot(xs,ys,'--g');hold on;
auc_arr(4)=auc;
s4 = sprintf('Gleason 3+4, AUC = %1.2f',auc);

load('43cf.mat');
plot(xs,ys,'r');hold on;
auc_arr(5)=auc;
s5 = sprintf('Gleason 4+3, AUC = %1.2f',auc);

load('44cf.mat');
plot(xs,ys,'b');hold on;
auc_arr(6)=auc;
s6 = sprintf('Gleason 4+4, AUC = %1.2f',auc);

load('45cf.mat');
plot(xs,ys,'k');hold on;
auc_arr(7)=auc;
s7 = sprintf('Gleason 4+5, AUC = %1.2f',auc);

load('53cf.mat');
plot(xs,ys,'g');hold on;
auc_arr(8)=auc;
s8 = sprintf('Gleason 5+3, AUC = %1.2f',auc);

load('54cf.mat');
plot(xs,ys,'-.r');hold on;
auc_arr(9)=auc;
s9 = sprintf('Gleason 5+4, AUC = %1.2f',auc);

load('55cf.mat');
plot(xs,ys,'-.b');hold on;
auc_arr(10)=auc;
s10 = sprintf('Gleason 5+5, AUC = %1.2f',auc);

load('BPcf.mat');
plot(xs,ys,':r');hold on;
auc_arr(11)=auc;
s11 = sprintf('BPH, AUC = %1.2f',auc);

load('NMcf.mat');
plot(xs,ys,':b');hold on;
auc_arr(12)=auc;
s12 = sprintf('Normal, AUC = %1.2f',auc);

load('HGcf.mat');
plot(xs,ys,':k');hold on;
auc_arr(12)=auc;
s13 = sprintf('HGPIN, AUC = %1.2f',auc);
h=xlabel('P_F(stroma)');
set(h,'FontSize',14);
h=ylabel('P_D(gland)');
set(h,'FontSize',14);
h=legend(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13);
set(h,'FontSize',14);
set(gca,'FontSize',12);
