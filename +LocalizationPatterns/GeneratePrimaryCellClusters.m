function GeneratePrimaryCellClusters(SettingsFilePath,BootstrapNumber)

% GeneratePrimaryCellClusters samples data for single cells according to 
%   parameters given by a settings file. It saves the clustering into a 
%   results file recording the bootstrap number.
%   
%   The setting file parameters are as follows:
%
%       - structSettings.PathToData ---> Path to file containing cell data.
%           e.g. C:\Path\To\Data\Data.mat
%       - structSettings.OutputFolder ---> Full path to output folder
%           e.g. C:\Output\folder\
%       - structSettings.MinClusterNumber ---> Minimal number of clusters.
%           e.g. 2
%       - structSettings.MaxClusterNumber ---> Maximal number of clusters.
%           e.g. 50
%       - structSettings.FirstRandomSampling  ---> Number of cells to be
%           sampled at random prior to calculation of data density. This is
%           used for data downsampling to avoid memory usage problems when
%           handling large number of cells. The number will depend on the 
%           memory availability of the current system. e.g. 70000
%       - structSettings.FeatureIndex  ---> Index of features in data set
%           to be used to generate the pattern classification. 
%           e.g. [1:5 7:10] 
%       - structSettings.NumberOfBins ---> Number of bins to divide the 
%           randomly sampled cells
%       - structSettings.CellSamplingsPerBin ---> Number of cells sampled
%           per bin
%
%
%   The output file:
%       - structResults.structSettings ---> The setting structure
%       - structResults.FullDataSize ---> The size of the input data
%       - structResults.SampleIndex1 ---> The sample idices for the first
%           random sampling
%       - structResults.SampleIndex2 ---> The sample indices for the second
%           random sampling
%       - structResults.matchix1 ---> Matching index of cells in first
%           sampling to cells of second samplings
%       - structResults.matchix2 ---> Matching index of cells in the second
%           sampling to cells in the first sampling
%       - structResults.z1 ---> Linkage results of first sampling
%       - structResults.z2 ---> Linkage results of second sampling
%       - structResults.IndexRep1 ---> Adjusted rand index of very 
%           clustering size for first sampling
%       - structResults.IndexRep2 ---> Adjusted rand index of very 
%           clustering size for first sampling
%
%
% Developed in University of Zurich, Institute of Molecular Life Sciences
% Copyright 2013.
%
% Authors:
%   Nico Battich
%   Thomas Stoeger
%   Lucas Pelkmans
%
% Website: https://www.pelkmanslab.org/
%




%--------------------------------------------------------------------------
% SETTINGS LOADING, DATA LOADING AND PARAMETERS CHECKS
fprintf('%s ::::::> Performing setting loadings and parameters checks\n',mfilename);

% load and evalute the settings file
[structSettings,CellData] = LocalizationPatterns.EvaluteSettingsFile(SettingsFilePath);
FullDataSize=size(CellData);

%--------------------------------------------------------------------------
% INITIALIZE RANDOM NUMBER GENERATOR (important for cluster computing)
rand('twister',sum(100*clock))
pause(BootstrapNumber/1000)
rand('twister',sum(100*clock))


%--------------------------------------------------------------------------
% INITIALIZE CLUSTERING ANALYSIS

% number of clusters
ClusterNumbers=[structSettings.MinClusterNumber:structSettings.MaxClusterNumber];
%structResults.ClusterNumbers=ClusterNumbers;

% features to be used 
FeatIx=structSettings.FeatureIndex;

% define number of cells to be sampled
BinNum=structSettings.NumberOfBins; 
SampleBin=structSettings.CellSamplingsPerBin; 

% generate two data samplings with replacement
tic
fprintf('%s ::::::> Subsampling data ::::::>',mfilename);
DownSample=structSettings.FirstRandomSampling;

% first sampling 
DownSamplingId=randperm(FullDataSize(1)); 
RandIx=DownSamplingId(1:DownSample);
DataRandomSample=CellData(RandIx,:);
EstimatedDensity=LocalizationPatterns.EstimateDataDensity(DataRandomSample(:,FeatIx),1.1,0.2);
SampleIndex1=LocalizationPatterns.SubSampleCells(RandIx,EstimatedDensity,BinNum,SampleBin);

% second sampling
DownSamplingId=randperm(FullDataSize(1));
RandIx=DownSamplingId(1:DownSample);
DataRandomSample=CellData(RandIx,:);
EstimatedDensity=LocalizationPatterns.EstimateDataDensity(DataRandomSample(:,FeatIx),1.1,0.2);
SampleIndex2=LocalizationPatterns.SubSampleCells(RandIx,EstimatedDensity,BinNum,SampleBin);

% clear data from memory
clear DataRandomSample DownSamplingId DownSample RandIx EstimatedDensity
toc

% get both random samplings
DataSamp1=CellData(SampleIndex1,FeatIx);
DataSamp2=CellData(SampleIndex2,FeatIx);

% clear data from memory
clear CellData


%--------------------------------------------------------------------------
% PERFORM CLUSTERING

% calculate cells closer neighboughrs, necesary for the adjusted rand index
% calculation
d1=pdist2(DataSamp1,DataSamp2);
[~,matchix1]=min(d1,[],2);
d2=pdist2(DataSamp2,DataSamp1);
[~,matchix2]=min(d2,[],2);
% clear data from memory
clear d1 d2

tic
fprintf('%s ::::::> Clustering data ::::::>',mfilename);

% distances (euclidean)
d1=pdist(nanzscore(DataSamp1));
d2=pdist(nanzscore(DataSamp2));

% perform linkage (ward, to get cluster of relativelly similar sizes)
try
z1=linkage(d1,'ward');
z2=linkage(d2,'ward');
catch
fprintf('')    
end
% get clusters id for all cluster numbers specified
clustid1=arrayfun(@(x) cluster(z1,'MaxClust',x),ClusterNumbers,'uniformoutput',false);
clustid2=arrayfun(@(x) cluster(z2,'MaxClust',x),ClusterNumbers,'uniformoutput',false);

% calculate the adjusted rand indexes for all cluster numbers, assuming
% that the second classification of a cell is that of its closest neighbour
matchid1=cellfun(@(x) x(matchix1),clustid2,'uniformoutput',false);
matchid2=cellfun(@(x) x(matchix2),clustid1,'uniformoutput',false);

IndexRep1=cellfun(@(a,b) LocalizationPatterns.adjrand(a,b),clustid1,matchid1);
IndexRep2=cellfun(@(a,b) LocalizationPatterns.adjrand(a,b),clustid2,matchid2);
toc


%--------------------------------------------------------------------------
% SAVE DATA

tic
fprintf('%s ::::::> Saving data...\n',mfilename);

% build results structure
structResults.structSettings=structSettings;
structResults.FullDataSize=FullDataSize;
structResults.SampleIndex1=SampleIndex1;
structResults.SampleIndex2=SampleIndex2;
structResults.matchix1=matchix1;
structResults.matchix2=matchix2;
structResults.z1=z1;
structResults.z2=z2;
structResults.IndexRep1=IndexRep1;
structResults.IndexRep2=IndexRep2;


strFileName=sprintf('ResultsPrimaryClustering_InterationNumber_%.4d.mat',BootstrapNumber);
%strfinalname=fullfile(structSettings.OutputFolder,strFileName);
if isunix
    strFileName=[structSettings.OutputFolder '/' strFileName];
else
    strFileName=[structSettings.OutputFolder '\' strFileName];    
end

save(strFileName,'structResults','-v7.3');

toc
