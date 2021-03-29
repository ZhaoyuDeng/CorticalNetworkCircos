% test function CircosDataOrganize.m
clc;clear;

% load raw data
load('RawDataCircos.mat');
% define variables
working_dir = pwd;
P_threshold = 0.005;
link_mode = 3;

% test function
CircosDataOrganize(working_dir,RawDataCircos,P_threshold,link_mode)


% run Circos command, if need run Matlab in Terminal
system('circos -conf CircosPlot.conf');