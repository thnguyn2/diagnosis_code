clc;
clear all;
close all;
%Get all the files name and of the glands and raw data
addpath(strcat(cd(cd('..')),'\support'));
datapath = 'H:\QPI data\';
[filenames,glandnames]=findFileName(datapath);
%Calculate the filteatapath);

%% Update the label of all classes
%updateLabelImages(filenames,glandnames,strcat(datapath,'label\'));

%% Calculate the phase distributiondist=updatePhaseDistribution(filenames,strcat(datapath,'phase_dist\'),strcat(datapath,'label\'));

%% response to odd & even symmetric filters
%updateEnergy(filenames,strcat(datapath,'texdir\'));%Calculate the response to the LM filter

%% update this histogram of texture directions- this is for HOG feature
%updateTextureDirectionDistribution(filenames,strcat(datapath,'label\'),strcat(datapath,'texdir\'),strcat(datapath,'text_dir_hist\'));

%% computing the texton map - This is for texton feature
%updateTexton(filenames,strcat(datapath,'label\'),strcat(datapath,'texdir\'));%Compute the texton map

%% update the histogram of texton distribution in a window
%updateTextonDistribution(filenames,strcat(datapath,'label\'),strcat(datapath,'texdir\'))

%% update affinity matrix for each data -  for FLD classifier (see old code)
%%%%%updateAffinityMatrix(filenames,strcat(datapath,'texdir\'));

%% compute the phase variance within a window - for phase variance matrix
%updatePhaseVariance(filenames,strcat(datapath,'label\'),strcat(datapath,'phase_dist\'))

%% update the weight - for normalized cut (see old code)
%%%%%updateWeight(filenames,strcat(datapath,'texdir\'),256);

%% calculate the filter response - for nuclei detection (see old code)
%%%%% updateFilterResForNucleiCalc(filenames,strcat(datapath,'ForNuclei\'));%Calculate the response to the LM filter

%% update the Hough transform - for nuclei dection (see old code)
%%%%%updateNuclei(filenames,strcat(datapath,'ForNuclei\'));

%% update the morphological information-From the groundtruth
%updateGlandSize(filenames,strcat(datapath,'label\'),strcat(datapath,'morpinfo\glandsize\'));

%% update the morphological information without groundtruth
%updateGlandSizeNoGT(filenames,strcat(datapath,'svm\texton_ws40_1024\'),strcat(datapath,'morpinfo\glandsize\'));
