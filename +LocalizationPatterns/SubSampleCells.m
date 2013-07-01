function SampledIndexes=SubSampleCells(RandIx,N,BinNum,SampleBin)

%BinNum=200;
BinBounds=min(N):(max(N)/(BinNum-1)):max(N);
[BinSize,BinId]=histc(N,BinBounds);
%SampleBin=40;

SampledIndexes=[];
for iBin=1:max(BinId)%length(unique(BinId))
    binix=find(BinId==iBin);
    if length(binix)<SampleBin
        SampledIndexes=[SampledIndexes;RandIx(binix)'];
    else
        tmpix=randperm(length(binix));
        SampledIndexes=[SampledIndexes;RandIx(binix(tmpix(1:SampleBin)))'];
    end
end
