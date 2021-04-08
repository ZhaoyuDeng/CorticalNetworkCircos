% test function CircosDataOrganize.m
clc;clear;

% load raw data
load('RawDataCircos.mat');
% define variables
workingDir = pwd;
LINK_MODE = 4;
% check existence of labels
if isfield(RawDataCircos,'HigherOrderNetworkLabel') && isfield(RawDataCircos,'ElementLabel')
    command = '';
else
    command = 'LT';
end

% test function
[fileFullPath1,fileFullPath2,fileFullPath3] = CircosDataOrganize(workingDir,RawDataCircos,LINK_MODE);
confFullPath = EditConf(workingDir,command);


% run Circos command, if need run Matlab in Terminal
system('circos -conf CircosPlot.conf');

%% codes organize new processed matrix
% load('RawDataCircos.mat');
% P_THRESHOLD = 0.005;
% % generate new pattern Matrix
% filMatCorr = RawDataCircos.Matrix;
% filMatCorr(RawDataCircos.P_Corrected > P_THRESHOLD) = 0;
% RawDataCircos.prosMatrix = filMatCorr;