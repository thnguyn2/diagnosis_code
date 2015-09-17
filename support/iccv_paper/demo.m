%% Update the label of all classesclc;
clc;
clear all;
close all;
%Get all the files name and of the glands and raw data
addpath(strcat(cd(cd('..')),'\support'));
datapath = 'H:\TMA_cores_and_diagnosis\';

%[filenames,glandnames]=findFileName(datapath);%This function is for the
%QPI_data folder
disp('Finding the list of all labels');
[filenames,glandnames,gradelist]=findFileNameFromROIs(datapath);
%save(strcat(datapath,'diag\fileinfo.mat'),'filenames','glandnames','gradelist');

%Calculate the filteatapath);


%disp('Updating labels for all images');
updateLabelImages(filenames,glandnames,strcat(datapath,'label\'));

%% Calculate the phase distributiondist=updatePhaseDistribution(filenames,strcat(datapath,'phase_dist\'),strcat(datapath,'label\'));

%% response to odd & even symmetric filters
%disp('Computing the filter response');
%updateEnergy(filenames,strcat(datapath,'texdir\'));%Calculate the response to the LM filter

%% [Non-texton--no need to run this]update this histogram of texture directions- this is for HOG feature
%updateTextureDirectionDistribution(filenames,strcat(datapath,'label\'),strcat(datapath,'texdir\'),strcat(datapath,'text_dir_hist\'));

%% computing the texton map - This is for texton feature
%updateTexton(filenames,strcat(datapath,'label\'),strcat(datapath,'texdir\'));%Compute the texton map

%% update the histogram of texton distribution in a window
updateTextonDistribution(filenames,strcat(datapath,'label\'),strcat(datapath,'texdir\'))
%
%% compute the phase variance within a window - for phase variance matrix
%updatePhaseVariance(filenames,strcat(datapath,'label\'),strcat(datapath,'phase_dist\'))

%% update the weight - for normalized cut (see old code)
%%%%%updateWeight(filenames,strcat(datapath,'texdir\'),256);

%% calculate the filter response - for nuclei detection (see old code)
%%%%% updateFilterResForNucleiCalc(filenames,strcat(datapath,'ForNuclei\'));%Calculate the response to the LM filter

%% update the Hough transform - for nuclei dection (see old code)
%%%%%updateNuclei(filenames,strcat(datapath,'ForNuclei\'));

%% update the morphological information-From the groundtruth
% updateGlandSize(filenames,strcat(datapath,'label\'),strcat(datapath,'diag\'),strcat(datapath,'texdir\'));%Compute the diagnosis information and save it

%% update the morphological information without groundtruth
%%%%%%updateGlandSizeNoGT(filenames,strcat(datapath,'svm\texton_ws40_1024\'),strcat(datapath,'morpinfo\glandsize\'));

%Compute the histogram information on all the image
%computeTextonHistForWholeCore(filenames,strcat(datapath,'label\'),strcat(datapath,'texdir\'));