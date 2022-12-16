%% Script to run microstate analysis. 
% note that the parameters given in this script have to be commmented out
% in the corresponding matlab files.
%% Preparation
clc, clear all, close all,
% Setup to get latex figures etc.
% run setup.m
preprocessing = 'preprocessed_220608_Filt_1_40';

%% Microstate analysis
tempfile = 'Literature';
Nmicrostates = [8:10];

%%
run microstate_individual_clustering_backfitting_rs_0807;

%% Entropy
date_today = '220609';
for i = 1:length(Nmicrostates)
    Nmicro = Nmicrostates(i);
    run entropy_individual_rs
end
