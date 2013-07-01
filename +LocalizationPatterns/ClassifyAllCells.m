function CellClassification=ClassifyAllCells(strPathToData,strPathToCentroidCluster,CentroidSamplingSize,BootstrapNumber)

% ClassifyAllCells classify single cells using the results of the centroid
% clustering from the function BuildCentroidCluster. 
%  
% Inputs are as follows:
%   - strPathToData ---> Path and fullname to cell data. e.g.
%       C:\Path\To\Data\Data.mat
%   - strPathToCentroidCluster ---> Path to the results from the
%       centroid clustering. e.g  C:\Path\To\Centroid\Clustering\Results.mat
%   - CentroidSamplingSize ---> Number of centroids to be sampled at each 
%       itteration of the classification. 
%       This is an optional parameter, default value is 10.
%   - BootstrapNumber ---> Number of bootstraps for classification of
%       cells. This is an optional parameter, default value is 1000.
%
%   Output:
%   - CellClassification ---> Classification weight or penetrance for 
%       each cell to each pattern
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
% CHECK INPUTS & LOAD DATA

if nargin<2
    error('Number of inputs not correct, the path to the cell data mut be given as well as the path to the results of the centroid clusterig.')
elseif nargin<3
    warning('Sampling size for centroids was not provided. Using default setting of 10 sampled centroids per itteration in classification.')
    CentroidSamplingSize=5;
    warning('Bootstrap size not provided. Using default setting of 1000 bootstraps.')
    BootstrapNumber=1000;
elseif nargin<4
    warning('Bootstrap size not provided. Using default setting of 1000 bootstraps.')
    BootstrapNumber=1000;
end


%--------------------------------------------------------------------------
% LOAD DATA
tic
fprintf('%s ::::::> Loading data ::::::>',mfilename);
% load data
load(strPathToData)
load(strPathToCentroidCluster)
toc
CellData=structDataToUse.CellData(:,structCentroidCluster.FeatureIndex);

%--------------------------------------------------------------------------
% INITIALIZE RANDOM NUMBER GENERATOR (important for cluster computing)
rand('twister',sum(100*clock))
pause(CentroidSamplingSize/1000)
rand('twister',sum(100*clock))


%--------------------------------------------------------------------------
% PERFORM CLASSIFICATION
tic
fprintf('%s ::::::> Classifying cells, pleas wait....\n',mfilename);

ClusterIds=[1:max(structCentroidCluster.ClusterId)]';
CellNumber=length(CellData);
CentroidNumber=length(structCentroidCluster.MeasuredCentroids);
PreCellClassification=nan(CellNumber,BootstrapNumber);

for i=1:BootstrapNumber
    randix=randperm(CentroidNumber);
    TemporaryClusterId=structCentroidCluster.ClusterId(randix(1:CentroidSamplingSize),1);
    d=pdist2(CellData,structCentroidCluster.MeasuredCentroids(randix(1:CentroidSamplingSize),:),'euclidean');
    [~,minix]=min(d,[],2);
    PreCellClassification(:,i)=TemporaryClusterId(minix);
end
fprintf('%s ::::::> ',mfilename);
toc

CellClassification=arrayfun(@(x) sum(PreCellClassification==x,2)./size(PreCellClassification,2),ClusterIds,'uniformoutput',false);
CellClassification=cat(2,CellClassification{:});
