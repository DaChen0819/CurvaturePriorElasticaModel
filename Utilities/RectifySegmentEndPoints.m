function [rectifiedPaths,rectifiedPoints]=RectifySegmentEndPoints(voronoiMap,taggedTrajectory,cost,options)
    
assert(nargin==4);


numOriens=options.numOriens;
taggedTrajectory=round(taggedTrajectory);
numTrajectories=round(max(voronoiMap(:)));
gridScale=options.gridScale;
origin=options.origin;
admPercent=options.admPercent;

options.euclideanScale=[1;1;0];
options.exportValues=1.0;
options.gridScale=gridScale;
options.origin=origin;
options.imageSize=size(taggedTrajectory);
options.model=options.model;
options.numOriens=numOriens;

[numRows,numCols]=size(taggedTrajectory);
input.model=options.model;
input.order=1.0;
input.seedRadius=0;
input.eps=0.1;
input.dims = [numCols;numRows;numOriens];
input.arrayOrdering='YXZ_ColumnMajor';
input.verbosity=0;
input.origin=options.origin;  % Physical origin
input.gridScale=options.gridScale; % Physical gridScale
input.exportValues=1.0;
input.exportRectifiedPoint=1.0;
input.euclideanScale=options.euclideanScale;
rectifiedPaths=cell(1,numTrajectories);
rectifiedPoints=[];

input.geodesicSolver='ODE';
input.geodesicStep=0.5;
input.geodesicCausalityTolerance=4;
input.geodesicTargetTolerance=6;
input.cost=cost;
input.model=options.model;
input.xi=options.xi;

for j=1:numTrajectories
    endPointsMap = bwmorph(taggedTrajectory==j,'endpoints');
    [mRows,mCols]=find(endPointsMap==1);
    mEndPoints=[mCols';mRows'];

    cEndPoints= RescaledCoords(mEndPoints,[0;0],[gridScale;gridScale],false);
    cSeeds=repmat(cEndPoints(:,1),[1 numOriens]);
    cOLSeeds=cat(1,cSeeds,(2*pi/numOriens).*(0:numOriens-1));% consider all the possible orientations.
    cTips=repmat(cEndPoints(:,2),[1 numOriens]);
    cOLTips=cat(1,cTips,(2*pi/numOriens).*(0:numOriens-1)); % consider all the possible orientations.

    options.mask=double(voronoiMap==j);
    trajectoryLength=TrajectoryLengthComputation(cOLSeeds,cOLTips,cost,options);

    input.physicalMask=double(voronoiMap==j);
    input.seeds=cOLSeeds;
    input.admEuclidLength=admPercent*trajectoryLength;
    
    output=HamiltonRectifySegmentEndPoints_CurvatureRegularization(input);
    rectifiedIndex=output.rectifiedIndex+1; % in C++ the index  is from 0 while in Matlab it is from 1.
    orien=(2*pi/numOriens)*(rectifiedIndex(3)-1);
    clear output;
    cRectifiedSeed=RescaledCoords(rectifiedIndex(1:2),origin,[gridScale;gridScale],false);
    cOLSeed=cat(1,cRectifiedSeed,orien+pi);
    input.seeds=cOLSeed;
    input.admEuclidLength=admPercent*trajectoryLength;
    output=HamiltonRectifySegmentEndPoints_CurvatureRegularization(input);

    rectifiedIndex=output.rectifiedIndex+1; %in C++ the index  is from 0 while in Matlab it is from 1.
    cGeoPaths_rectified=mat2cell(output.geodesicPoints_rectified(1:3,:),3,output.geodesicLengths_rectified);
    rectifiedPaths{j}=cGeoPaths_rectified{1};
    rectifiedPoints=cat(2,rectifiedPoints,rectifiedIndex);
end

end

function trajectoryLength=TrajectoryLengthComputation(seeds,tips,cost,options)
    numRows=options.imageSize(1);      
    numCols=options.imageSize(2);
    numOriens=options.numOriens;
    input.euclideanScale=options.euclideanScale;
    input.dims=[numCols;numRows;numOriens];
    input.gridScale=options.gridScale;
    input.origin=[0.0;0.0];
    input.exportValues=options.exportValues; 

    input.factoringMethod='None';
    input.verbosity=0;
    input.arrayOrdering='YXZ_ColumnMajor';
    input.order=1;
    input.exportGeodesicFlow=0.0;
    input.exportEuclideanLengths=1.0;
    input.model=options.model;
    input.cost=cost;
    input.seeds=seeds;
    input.stopWhenAnyAccepted=tips;

    input.activeRegion=repmat(options.mask,[1 1 numOriens]);
    output_GeoDistance=HamiltonGeodesicDistanceComputation(input);
    mFirstRearchedTip=output_GeoDistance.firstAcceptedTip+1.0;
    trajectoryLength=output_GeoDistance.euclideanLengths...
        (mFirstRearchedTip(2),mFirstRearchedTip(1),mFirstRearchedTip(3));
end





