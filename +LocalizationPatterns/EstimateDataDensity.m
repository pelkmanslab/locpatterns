function DataDensity=EstimateDataDensity(Data,CorrectionFactor,divFac)

% divFac must be between 0 and 1, it depends on the amount of memory
% available and speed required
% CorrectionFactor is used to correct the radious for calculation of number
% of neiboughrs recomended between 0.9 and 1.1

% initialize parameters 
obnum=size(Data,1);
DataDensity=nan(obnum,1);

% do density estimation
randix=randperm(obnum);
d=pdist2(Data,Data(randix(1:round(obnum.*divFac)),:));
d(d==0)=nan;
Radious=quantile(min(d,[],2),0.999);
Radious=CorrectionFactor*Radious;
d=d<Radious;
n=sum(d,2);
DataDensity=n;
%figure;hist(DataDensity,100)
end