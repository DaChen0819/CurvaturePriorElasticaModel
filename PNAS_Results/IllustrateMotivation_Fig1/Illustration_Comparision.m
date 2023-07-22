clear variables; close all; clc;

load('rawImg.mat');
normalizedImg=im2double(rawImg(:,:,1));
[numRows,numCols]=size(normalizedImg);

%% original setting.
numOriens=72;
gridScale=2*pi/numOriens;
origin=[0;0];
maximalRadius_Rectified=2;
maximalRadius_Tube=3;
BrightObject=1;
alpha=4;
minBranchLength=0;
beta_Data=6.5;
beta_xi=3;
mTipPt=[274;26];
mTipTangentPt=[275;26];
mSrcPt=[40;15];
mSrcTangentPt=[40;12];

%% points.
theta_Src=TransferTangentToOrientation(mSrcTangentPt-mSrcPt);
theta_Tip=TransferTangentToOrientation(mTipTangentPt-mTipPt);
cRescaledSeeds=RescaledCoords(mSrcPt,origin,[gridScale;gridScale],false);
input_PriorGeos.seeds=cat(2,[cRescaledSeeds;theta_Src],[cRescaledSeeds;theta_Src+pi]);
cRescaledTips=RescaledCoords(mTipPt,origin,[gridScale;gridScale],false);
input_PriorGeos.stopWhenAnyAccepted=cat(2,[cRescaledTips;theta_Tip],[cRescaledTips;theta_Tip+pi]);

%% compute orientation scores
try
input_OOF.InputImage=double(normalizedImg);
input_OOF.ComputeOrientationScore=1; 
input_OOF.ComputeMultiScaleTubularScore=1; 
input_OOF.MultiScaleEigenAnalysis=1;
input_OOF.BrightObject=BrightObject; 
input_OOF.FixedSigma=1;
input_OOF.MinSigma=1;
input_OOF.MaxSigma=4; 
input_OOF.NBOfScales=2*(input_OOF.MaxSigma-input_OOF.MinSigma)+1;
input_OOF.NBOfOrientations=numOriens;
input_OOF.InitialOrien=[0.0;-1.0];
OOF=OptimallyOrientedFluxFilter(input_OOF);
orienScore=OOF.OrienScore;
normalizedOS=orienScore./(max(orienScore(:)));
clear orienScore; clear OOF;

%% compute the centerlines of vessel segments.
thresValue=0.4;
binarySeg=normalizedImg>thresValue;
binarySeg=bwmorph(binarySeg,'open');
binarySeg=imfill(binarySeg,'holes');
skeletonMap = bwskel(binarySeg,'MinBranchLength',minBranchLength);
branchPtMap=bwmorph(skeletonMap,'branchpoints');
newSkeletonMap=RemoveNeighPointsFromEnlargedNeighbourhood(double(skeletonMap),double(branchPtMap));
seperateSkeletonMap=RemoveSmallSegments(newSkeletonMap,5);
seperateSkeletonMap=RemoveCircularSegments(seperateSkeletonMap);
%% compute the curvature of these vessel centerlines.
options.eps=0.1;
options.numOriens=numOriens;
options.euclideanScale=[1.0;1.0;0.0];
options.model='Elastica2';
options.gridScale=gridScale;
options.origin=origin;
options.admPercent=1.0;
options.xi=1;
options.maximalRadius=maximalRadius_Rectified;
options.imageSize=[numRows;numCols];
[voronoiMap,taggedTrajectory]=ConstructVoronoiTessellation(seperateSkeletonMap,options.maximalRadius);
[cPaths_Rectified,~]=RectifySegmentEndPoints(voronoiMap,taggedTrajectory,exp(-alpha.*normalizedOS),options);
options.maximalRadius=maximalRadius_Tube;
tubularNeighRegion=TubularNeighbourhood_OpenCurve(cPaths_Rectified,options);

%% geodesic models.
input_PriorGeos.order=1.0;
input_PriorGeos.seedRadius=0;
input_PriorGeos.eps=0.1;
input_PriorGeos.dims = [numCols;numRows;numOriens];
input_PriorGeos.arrayOrdering='YXZ_ColumnMajor';
input_PriorGeos.verbosity=0;
input_PriorGeos.origin=origin;  % Physical origin
input_PriorGeos.gridScale=gridScale; % Physical gridScale
input_PriorGeos.exportValues=1; % distance table, of size [n,n,numberOfDirections]
input_PriorGeos.model='ElasticaExt2'; %Alternatively : 'ReedsSheppExt2','ReedsSheppForwardExt2', 'ElasticaExt2', 'DubinsExt2'
input_PriorGeos.geodesicSolver='ODE';
input_PriorGeos.geodesicStep=0.5;
input_PriorGeos.geodesicCausalityTolerance=4;
input_PriorGeos.geodesicTargetTolerance=6;

options.smoothWidth=5;
[~,extendedVelocity]=ConstructCurvaturePriorMap(cPaths_Rectified,options,double(tubularNeighRegion>0.5));
input_PriorGeos.xi =beta_xi;
input_PriorGeos.kappa=extendedVelocity;
input_PriorGeos.cost = exp(-beta_Data.*normalizedOS); 
output_Prior=HamiltonCurvaturePriorMinimalPaths(input_PriorGeos);
cGeodesics_Prior = mat2cell(output_Prior.geodesicPoints,3,output_Prior.geodesicLengths); 

%% classical elastica model.
input_PriorGeos.kappa=zeros(numRows,numCols,numOriens);
output_Elastica=HamiltonCurvaturePriorMinimalPaths(input_PriorGeos);
cGeodesics_Elastica = mat2cell(output_Elastica.geodesicPoints,3,output_Elastica.geodesicLengths); 

mRescaledGeo_Elastica=RescaledCoords(cGeodesics_Elastica{1}(1:2,:),origin,[gridScale;gridScale]);
mRescaledGeo_Prior=RescaledCoords(cGeodesics_Prior{1}(1:2,:),origin,[gridScale;gridScale]);
catch
    load(strcat(dataDir,'mRescaledGeo_Elastica'));
    load(strcat(dataDir,'mRescaledGeo_Prior'));
end


figure(1);imagesc(rawImg); axis off;axis image;colormap gray;hold on;
plot(mRescaledGeo_Prior(1,:),mRescaledGeo_Prior(2,:),'Color','r','LineWidth',3);
plot(mSrcPt(1),mSrcPt(2),'o','MarkerEdgeColor','r','MarkerFaceColor','r','markersize',10);
plot(mTipPt(1),mTipPt(2),'o','MarkerEdgeColor','b','MarkerFaceColor','b','markersize',10);


figure(2);imagesc(rawImg);axis off; axis image;hold on;colormap gray;
plot(mRescaledGeo_Elastica(1,:),mRescaledGeo_Elastica(2,:),'r','LineWidth',3);
plot(mSrcPt(1),mSrcPt(2),'o','MarkerEdgeColor','r','MarkerFaceColor','r','markersize',10);
plot(mTipPt(1),mTipPt(2),'o','MarkerEdgeColor','b','MarkerFaceColor','b','markersize',10);



