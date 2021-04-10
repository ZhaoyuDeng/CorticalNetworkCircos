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
[filePathBand,filePathLabel,filePathLink] = CircosDataOrganize(workingDir,RawDataCircos,LINK_MODE);
filePathConf = EditConf(workingDir,command);


% run Circos command, if need run Matlab in Terminal
system('circos -conf CircosPlot.conf');

