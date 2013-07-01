function NameOfResultFile=BuildCentroidCluster(strPathToData,strPathToPrimaryClusterResults,ClusterNumber)

% BuildCentroidCluster takes the results from the primary cell clustering
% obtained by the GeneratePrimaryCellClusters funtion of the
% LocalizationPatterns module. It saves the centroid clustering output that
% then can be used to classify cells with the ClassifyAllCells function.
%
% Inputs are as follows:
%   - strPathToData ---> Path and fullname to cell data. e.g.
%       C:\Path\To\Data\Data.mat
%   - strPathToPrimaryClusterResults ---> Path to the results from the
%       primary clustering are located. The folder will be used as output
%       folder for the result file generated.
%   - ClusterNumber ---> Number of clusters cenroids should be divided
%       into. This is an optional parameter, default value is 5.
% 
% Output file:
%   - structCentroidCluster.z ---> Cluster final linkage
%   - structCentroidCluster.ClusterId ---> Cluster Id for each centroid
%   - structCentroidCluster.MeasuredCentroids ---> Measured centroid in 
%       feature space
%   - structCentroidCluster.FeatureIndex ---> Index of features used
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
    error('Number of inputs not correct, the path to the cell data mut be given as well as the path to the results of the primary cell clustering.')
elseif nargin<3
    warning('Final number of clusters not provided. Using default setting of 5 clusters.')
    ClusterNumber=5;
end

% load data 
load(strPathToData);

% get result files from primary clustering
FilesInDir=dir(strPathToPrimaryClusterResults);
FilesInDir=cat(1,{FilesInDir(:).name});

% the name of the primary results file
TestingExpression='ResultsPrimaryClustering_InterationNumber_.*.mat';
FilesFilter=cellfun(@(x) ~isempty(regexp(x,TestingExpression)),FilesInDir);
FilesInDir(~FilesFilter)=[];
NumberOfFiles=length(FilesInDir);

%--------------------------------------------------------------------------
% CALCULATE CENTROIDS

MeasuredCentroids=[];

for i=1:NumberOfFiles
    
    fprintf('%s: Dealing with file %.10d out of %.10d.\n',mfilename,i,NumberOfFiles);
    
    % load primary clustering results  
    load(fullfile(strPathToPrimaryClusterResults,FilesInDir{i}));
    
    % get the temporary cell data
    TemporaryData=structDataToUse.CellData(structResults.SampleIndex1,structResults.structSettings.FeatureIndex);
    
    % recluster linkage results
    TemporaryClusterIndex=cluster(structResults.z1,'MaxClust',ClusterNumber);
    
    % divide data according to cluster index
    GroupedMeanClusters=arrayfun(@(x) nanmean(TemporaryData(TemporaryClusterIndex==x,:),1),[1:ClusterNumber]','uniformoutput',false);
    
    % if number of cells are small or randon numbber generator is not 
    % initialized in a good way there can be repetition of samplings 
    if i>1
        TestRandomSampling=~ismember(GroupedMeanClusters{1},MeasuredCentroids,'rows');
    else
        TestRandomSampling=true;
    end
    
    if TestRandomSampling        
        MeasuredCentroids=[MeasuredCentroids;cat(1,GroupedMeanClusters{:})];
    end

    % repeat for second sampling
    % get the temporary cell data
    TemporaryData=structDataToUse.CellData(structResults.SampleIndex2,structResults.structSettings.FeatureIndex);
    
    % recluster linkage results
    TemporaryClusterIndex=cluster(structResults.z2,'MaxClust',ClusterNumber);
    
    % divide data according to cluster index
    GroupedMeanClusters=arrayfun(@(x) nanmean(TemporaryData(TemporaryClusterIndex==x,:),1),[1:ClusterNumber]','uniformoutput',false);
    
    % if number of cells are small or randon numbber generator is not 
    % initialized in a good way there can be repetition of samplings 
    if i>1
        TestRandomSampling=~ismember(GroupedMeanClusters{1},MeasuredCentroids,'rows');
    else
        TestRandomSampling=true;
    end
    
    if TestRandomSampling        
        MeasuredCentroids=[MeasuredCentroids;cat(1,GroupedMeanClusters{:})];
    end
    
end


%--------------------------------------------------------------------------
% CLUSTERING OF CENTROIDS

d=pdist(nanzscore(MeasuredCentroids));
z=linkage(d,'ward');
cid=cluster(z,'MaxClust',ClusterNumber);

%--------------------------------------------------------------------------
% SAVE DATA

structCentroidCluster.z=z;
structCentroidCluster.ClusterId=cid;
structCentroidCluster.MeasuredCentroids=MeasuredCentroids;
structCentroidCluster.FeatureIndex=structResults.structSettings.FeatureIndex;

NameOfResultFile=sprintf('ResultsCentroidCluster_ClusterNumber_%.4d.mat',ClusterNumber);

if isunix
    strFileName=[strPathToPrimaryClusterResults '/' NameOfResultFile];
else
    strFileName=[strPathToPrimaryClusterResults '\' NameOfResultFile];    
end

save(strFileName,'structCentroidCluster','-v7.3')

