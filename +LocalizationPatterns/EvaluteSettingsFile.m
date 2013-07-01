function [structSettings,CellData] = EvaluteSettingsFile(SettingsFilePath)

% EvaluteSettingsFile evaluates the setting file for pymary cell
% clustering.
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

fid=fopen(SettingsFilePath);
% loop over each line
fprintf('%s: Parsing setting file %s.\n',mfilename,SettingsFilePath);
while 1
    tline = fgetl(fid);
    if tline==-1,   break,   end
    
    %%% hehe, this is higly insecure! ALMOST arbitrary code
    %%% execution, woot! :D
    tline = strtrim(tline);
    
    if strncmpi(tline, 'structSettings.',15)
        
        % A setting was found. The default values will be changed
        strFieldName = regexpi(tline,'^structSettings.(\w{1,})','Tokens');
        strFieldName = char(strFieldName{:});
        
        if strncmp(strFieldName, 'PathToData',10)...
                || strncmp(strFieldName, 'ObjectName',10)...
                || strncmp(strFieldName, 'OutputFolder',12)...
                || strncmp(strFieldName, 'MinClusterNumber',16)...
                || strncmp(strFieldName, 'MaxClusterNumber',16)...
                || strncmp(strFieldName, 'FirstRandomSampling',19)...
                || strncmp(strFieldName, 'FeatureIndex',12)...
                || strncmp(strFieldName, 'NumberOfBins',12)...
                || strncmp(strFieldName, 'CellSamplingsPerBin',19)
            
            eval(tline);
            
        else
            error('%s: Typing error at %s',mfilename,tline)
        end
    end
end


% make sure all fields of the settings file exist, else provide default
% values

% check path to data
if ~isfield(structSettings,'PathToData')
    error('No data provided. Make sure that the settings file contains\na "PathToData" field.')
else
    %load data structDataToUse.fMetaData=fMetaData;
    fprintf('%s ::::::> Loading data\n',mfilename);
    load(structSettings.PathToData);
    CellData=structDataToUse.CellData; 
    FullDataSize=size(CellData);
    clear structDataToUse
end

% check path to output folder
if ~isfield(structSettings,'OutputFolder')
    error('No output folder provided. Make sure that the settings file contains\na "OutputFolder" field.')
elseif ~isdir(structSettings.OutputFolder)
    fprintf('%s ::::::> Creating output folder\n',mfilename);
    mkdir(structSettings.OutputFolder);
end



% check minimal number of clusters
if ~isfield(structSettings,'MinClusterNumber')
    warning('No minimal number of clusters provided.\n---> MinClusterNumber set to 2n\')
    structSettings.MinClusterNumber=2;
elseif structSettings.MinClusterNumber<2
    warning('Minimal number of clusters provided is less than 2.\n--->MinClusterNumber set to 2n\')
    structSettings.MinClusterNumber=2;
end

% check maximal number of clusters
if ~isfield(structSettings,'MaxClusterNumber')
    warning('No maximal number of clusters provided.\n---> MaxClusterNumber set to 10n\')
    structSettings.MaxClusterNumber=10;
end

% check number of cells to be used for down sampling
if ~isfield(structSettings,'FirstRandomSampling')
    warning('No random sampling number provided.\nNote that this is a very important parameter as\na low number might introduce classification errors\n---> Gessing FirstRandomSamplingn\')
    if FullDataSize(1)>20000 && FullDataSize(1)<140000
        structSettings.FirstRandomSampling=round(FullDataSize(1)/2);
    elseif FullDataSize(1)<20000
        structSettings.FirstRandomSampling=FullDataSize(1);
    else
        structSettings.FirstRandomSampling=70000;
    end
end

% check which features should be used
if ~isfield(structSettings,'FeatureIndex')
    warning('No features to be used were spesified.\n---> All features will be used\n')
    structSettings.FeatureIndex=1:FullDataSize(2);
end

% check number of bins for the data to be divided
if ~isfield(structSettings,'NumberOfBins')
    warning('No number of bins for data division provided.\n--->NumberOfBins set to 50n\')
    structSettings.NumberOfBins=50;
end

% check number of cells to be sampled per bin
if ~isfield(structSettings,'CellSamplingsPerBin')
    warning('No number of cells sampled per bin was provided.\n--->CellSamplingsPerBin set to 100n\')
    structSettings.CellSamplingsPerBin=100;
end

% test min and max cluster numbers
if structSettings.MinClusterNumber>structSettings.MaxClusterNumber
    error('Please make the minimum number of clusters smaller than the maximum number of clusters.')
end
