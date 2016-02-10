%Display all ROC curves for different features
clc;
clear all;
close all;

%Draw the ROC curve for mean and std of the phase
load('avgstd.mat');
nradius = 7;
figure(1);
for radidx = nradius:nradius
    figure(1);
    plot(xs_std{radidx},ys_std{radidx},'k','linewidth',2);
    drawnow;
    hold on
    plot(xs_mean{radidx},ys_mean{radidx},'g','linewidth',2);
    hold on
    drawnow;
end
   
load('oehod.mat');
nradius = 5;
for radidx = nradius:nradius
    figure(1);
    plot(xs_oe{radidx},ys_oe{radidx},'b','linewidth',2);
    hold on;
    drawnow;
end
figure(1);
plot(xs_hod,ys_hod,'m','linewidth',2);
hold on;

load('textonroc.mat');
figure(1);
plot(xs_svm{2},ys_svm{2},'r','linewidth',2);
      
legend('PD','PM','OE','HOD','Texton');
xlabel('P_F(stroma)');
ylabel('P_D(gland)');

