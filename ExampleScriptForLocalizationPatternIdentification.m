% This script file contains an example of how to de the localization
% pattern analysis presented in Battich et al., 2013.
% 
% The example data provided consists of per cell localization features for 
% five genes. Example data was normailized by z-scoring. The rows in the 
% dataset represent cells, and colums the 32 per cell localization features
% discussed by Battich et al., 2013. 
%
% 
% TO RUN THIS SCRIPT YOU MUST:
%        1. Change the working directory to the directory where this script
%           and the LocalizationPatterns module are located
%        2. Change in the example settings file the absolute path to the 
%           example data
%        3. Define in the example settings file the output folder
%        4. Change in this script the absolute path of the settings file
%        5. Run the script
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
%


%--------------------------------------------------------------------------
% DEFINE EXAMPLE DATASET & OUTPUT FOLDERS

% define and load setting file for reference

% Absolute path to settimgs file
SettingsFilePath='C:\Path\To\Settings\File.txt';
settings=LocalizationPatterns.EvaluteSettingsFile(SettingsFilePath);


%--------------------------------------------------------------------------
% GENERATE PRIMARY CELL CLUSTERS

% primary celll clusters are generaed according to the settings file.
% look at help LocalizationPatterns.GeneratePrimaryCellClusters for more
% information.
ItterationNumber=30;
for i=1:ItterationNumber;
    LocalizationPatterns.GeneratePrimaryCellClusters(SettingsFilePath,i);
end

%--------------------------------------------------------------------------
% BUILD CENTROID CLUSTER

% define the final number of clusters. Note: Battich et al., used adjusted
% randidexes to define this number. You can find calculation of the
% adjusted rand indexes in the output file of
% LocalizationPatterns.GeneratePrimaryCellClusters
NumberOfFinalClusters=5;

% generate the centroids cluster
NameOfResultFile=LocalizationPatterns.BuildCentroidCluster(settings.PathToData,settings.OutputFolder,NumberOfFinalClusters);

%--------------------------------------------------------------------------
% CLASSIFY ALLL CELLS

% number of centroids to sampe per itteration 
CentroidSamplingSize=7;

% number of itterations
SecondItterationNumber=1500;

% define the file containing the centroid clustering. 
if isunix
    CentroisClusteringFile=[settings.OutputFolder '/' NameOfResultFile];  
else
    CentroisClusteringFile=[settings.OutputFolder '\' NameOfResultFile];  
end

% perform cell classification
CellClassification=LocalizationPatterns.ClassifyAllCells(settings.PathToData,CentroisClusteringFile,CentroidSamplingSize,SecondItterationNumber);

%--------------------------------------------------------------------------
% PLOT EXAMPLE RESULTS

% Note: all gene names correspond to the example data provided with this
% example script

load(settings.PathToData);
genes=unique(structDataToUse.Entrez)';
GeneNames={'LIPA','COX1','ND2','ND3','TFRC','CXCR4'};
x=0:0.04:1;
histograms=arrayfun(@(a) histc(CellClassification(structDataToUse.Entrez==a,:),x),genes,'uniformoutput',false);

colors={'r','g','b','k','c','m'};
pattern=1;
%figure;
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1500, 895]);
subplot(2,3,1)
for i=1:length(histograms)
plot(x,histograms{i}(:,pattern)./sum(histograms{i}(:,pattern)),colors{i})
hold on
end
title('patter 1')
ylabel('fraction of cells')
xlabel('classification penetrance')
h = legend(GeneNames,6);
set(h,'Interpreter','none')
pattern=2;
subplot(2,3,2)
for i=1:length(histograms)
plot(x,histograms{i}(:,pattern)./sum(histograms{i}(:,pattern)),colors{i})
hold on
end
title('patter 2')
ylabel('fraction of cells')
xlabel('classification penetrance')
pattern=3;
subplot(2,3,3)
for i=1:length(histograms)
plot(x,histograms{i}(:,pattern)./sum(histograms{i}(:,pattern)),colors{i})
hold on
end
title('patter 3')
ylabel('fraction of cells')
xlabel('classification penetrance')
pattern=4;
subplot(2,3,4)
for i=1:length(histograms)
plot(x,histograms{i}(:,pattern)./sum(histograms{i}(:,pattern)),colors{i})
hold on
end
title('patter 4')
ylabel('fraction of cells')
xlabel('classification penetrance')
pattern=5;
subplot(2,3,5)
for i=1:length(histograms)
plot(x,histograms{i}(:,pattern)./sum(histograms{i}(:,pattern)),colors{i})
hold on
end
title('patter 5')
ylabel('fraction of cells')
xlabel('classification penetrance')

% END OF PLOT & FILE
%--------------------------------------------------------------------------
