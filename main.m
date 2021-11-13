clear;
clc;
close all;
%%
addpath('./biodata');
addpath('./tool');
% Load the multi-omics datasets
importdata('./biodata/GBM/GLIO_Gene_Expression.txt');
X1=ans.data;
importdata('./biodata/GBM/GLIO_Methy_Expression.txt');
X2=ans.data;
importdata('./biodata/GBM/GLIO_Mirna_Expression.txt');
X3=ans.data;
importdata('./biodata/GBM/GLIO_Survival.txt');
SUR=ans.data;

X1=normalize(X1);
X2=normalize(X2);
X3=normalize(X3);
alldata{1}=X1;
alldata{2}=X2;
alldata{3}=X3;
DATA={alldata{1},alldata{2},alldata{3}};
% Run our proposed method
nclass = 3;
opts.clusternum=nclass;
opts.beta=4;
gamma=6;
ksk=ConstructA_3order_NaN(DATA,gamma,nclass);% Construct the high order similarity matrices 
[S,w] = MVMLV(ksk,opts);% Calculate the fusion similarity matrix
group = SpectralClustering(S,nclass);% Extract group information using spectral clustering method
[p,fh,stats]=MatSurv(SUR(:,1), SUR(:,2),cellstr(int2str(group)),'NoPlot',false);% Draw survival curves 
