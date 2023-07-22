function [voronoiMap,taggedTrajectory]=ConstructVoronoiTessellation(unLabeledSeedMap,maxRadius)


if nargin==1
    maxRadius=10;
end

[numRows,numCols,~]=size(unLabeledSeedMap);
unLabeledSeedMap=double(unLabeledSeedMap);
unLabeledSeedMap=unLabeledSeedMap./max(unLabeledSeedMap(:));

[linearizedSeeds,seedFlags,taggedTrajectory]=LinearizeSeeds(unLabeledSeedMap);
[subY,subX]=ind2sub([numRows,numCols],linearizedSeeds);

mSeeds=[subX;subY];
cSeeds=mSeeds-1.0;


% Defining the domain
input.order=1.0;
input.seedRadius=0;
input.eps=0.1;
input.dims = [numCols;numRows];
input.arrayOrdering='YXZ_ColumnMajor';
input.verbosity=0;
input.origin=[0;0];  
input.gridScale=1.0; 
input.stopAtDistance=maxRadius;
input.voronoiStoppingCriterion='None';
input.exportValues=1.0;
input.exportVoronoiFlags=1.0;
input.exportGeodesicFlow=0.0;
input.model='Isotropic2';
input.cost=1.0;
input.seeds=cSeeds;
input.seedFlags=seedFlags';
input.exportEuclideanLengths=0.0;

output=HamiltonGeodesicDistanceComputation(input);
voronoiMap=output.voronoiFlags;





end

function [linearizedSeedIdx,seedFlags,trajectoryLabelMap]=LinearizeSeeds(seedImage)
[numRows,numCols]=size(seedImage);
trajectoryLabelMap=zeros(numRows,numCols);

binaryImg=seedImage>0.5;
CC = bwconncomp(binaryImg);
linearizedSeedIdx=[];
seedFlags=[];
for j=1:CC.NumObjects
    PixelIdx=CC.PixelIdxList{j};
    PixelIdx=PixelIdx';
    linearizedSeedIdx=cat(2,linearizedSeedIdx,PixelIdx);
    trajectoryLabelMap(PixelIdx)=j;
    flags=j.*double(PixelIdx>0);
    seedFlags=cat(2,seedFlags,flags);
end


end

