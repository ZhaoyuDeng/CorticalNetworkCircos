% test function CircosDataOrganize.m
clc;clear;

% load raw data
load('RawDataCircos.mat');
% define variables
workingDir = pwd;
P_THRESHOLD = 0.005;
LINK_MODE = 4;

% test function
CircosDataOrganize(workingDir,RawDataCircos,P_THRESHOLD,LINK_MODE)


% run Circos command, if need run Matlab in Terminal
system('circos -conf CircosPlot.conf');