clear variables; clc; close all;

numGrids=191;
numOriens=72;
posSrc=1500; 
posEnd=3850;

% Defining the domain.
gridScale=2*pi/numOriens;
origin=[0.0;0.0];
mSpiral=MakeSpiralLine(3,6,4*pi+0.5*pi);
cRefPath=RescaledCoords(mSpiral(1:2,:),origin,[gridScale;gridScale],false);
cRefPath=cat(1,cRefPath,mSpiral(3,:));

figure(1);imagesc(ones(191,191),[0,1]);colormap gray;axis image;hold on;
plot(mSpiral(1,posSrc:posEnd),mSpiral(2,posSrc:posEnd),'r');

input.order=1.0;
input.seedRadius=0;
input.eps=0.1;
input.dims = [numGrids;numGrids;numOriens];
input.arrayOrdering='YXZ_ColumnMajor';
input.verbosity=0;
input.origin=origin; 
input.gridScale=gridScale; 
input.cost = 1; 

try
mSeedPos=RescaledCoords(cRefPath(1:2,posEnd),origin,[gridScale;gridScale]);
mSrcTangentPt=RescaledCoords(cRefPath(1:2,posEnd-1),origin,[gridScale;gridScale]);
seedOrien=TransferTangentToOrientation(mSrcTangentPt-mSeedPos);
input.seeds=cRefPath(1:3,posEnd);

mTipTangentPt=RescaledCoords(cRefPath(1:2,posSrc-1),origin,[gridScale;gridScale]);
mTipPos=RescaledCoords(cRefPath(1:2,posSrc),origin,[gridScale;gridScale]);
tipOrien=TransferTangentToOrientation(mTipTangentPt-mTipPos);
input.stopWhenAnyAccepted = cRefPath(1:3,posSrc); 

input.exportValues=1; % distance table, of size [n,n,numberOfDirections]
input.geodesicSolver='ODE';
input.geodesicStep=0.5;
input.geodesicCausalityTolerance=4;
input.geodesicTargetTolerance=6;

options.numOriens=numOriens;
options.gridScale=gridScale;
options.origin=origin;
options.imageSize=[numGrids;numGrids];
options.maximalRadius=1;
options.smoothWidth=5;
tubularNeighRegion=TubularNeighbourhood_OpenCurve({cRefPath},options);

input.xi = 10; 
[~,extendedVelocity]=ConstructCurvaturePriorMap({cRefPath},options,double(tubularNeighRegion>0.5));
input.model='ElasticaExt2'; 
input.kappa=extendedVelocity;
output_Prior=HamiltonCurvaturePriorMinimalPaths(input);
cGeodesics_Prior = mat2cell(output_Prior.geodesicPoints,3,output_Prior.geodesicLengths); 

catch
load('data/cGeodesics_Prior.mat');
end



mRescaledRefPath=RescaledCoords(cGeodesics_Prior{1}(1:2,:),input.origin,[input.gridScale;input.gridScale]);

l=22;
figure(2);imagesc(ones(numGrids,numGrids,3)); axis image;hold on;grid off;
mRescaledGeodesic1=RescaledCoords(cRefPath(1:2,posSrc:posEnd),origin,[gridScale;gridScale]);
plot(mRescaledGeodesic1(1,1:2:end),mRescaledGeodesic1(2,1:2:end),'r','LineWidth',2);
h5=func_DrawArrow(mSeedPos,mSeedPos+l*[cos(seedOrien);sin(seedOrien)],25,'BaseAngle',30,'Width',3); 
set(h5,'Facecolor','r','EdgeColor','r');
h6=func_DrawArrow(mTipPos,mTipPos+l*[cos(tipOrien);sin(tipOrien)],25,'BaseAngle',30,'Width',3);
set(h6,'Facecolor','b','EdgeColor','b');
plot(mSeedPos(1),mSeedPos(2),'o','MarkerEdgeColor','r','MarkerFaceColor','r','Markersize',10);
plot(mTipPos(1),mTipPos(2),'o','MarkerEdgeColor','b','MarkerFaceColor','b','Markersize',10);
xticklabels('');yticklabels('');
seedPosShifted=mSeedPos+[-8;-8];
tipShifted=mTipPos+[-12;1];
text(seedPosShifted(1),seedPosShifted(2),'$s$','FontSize',22,'Interpreter','latex');
text(tipShifted(1),tipShifted(2),'$y$','FontSize',22,'Interpreter','latex');
alpha(h5,0.3);


figure(3);imagesc(ones(numGrids,numGrids,3)); axis image;hold on;grid off;
plot(mRescaledRefPath(1,:),mRescaledRefPath(2,:),  'b','LineWidth',2);
plot(mRescaledGeodesic1(1,:),mRescaledGeodesic1(2,:),'--r','LineWidth',1);
xticklabels('');yticklabels('');
