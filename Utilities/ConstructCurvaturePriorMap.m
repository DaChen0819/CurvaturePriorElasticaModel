function [phyCurvatureMaps,extendedCurvatureVelocity]=ConstructCurvaturePriorMap(cOLPaths,options,mask)

numOriens=options.numOriens;
gridScale=options.gridScale;
origin=options.origin;
numRows=options.imageSize(1);
numCols=options.imageSize(2);

if nargin==2
    mask=ones(numRows,numCols);
end

curvature_phi=zeros(numRows,numCols);
orien_phi=zeros(numRows,numCols);
scoreMap_c=zeros(numRows,numCols);
scoreMap_o=zeros(numRows,numCols);
smoothWidth=options.smoothWidth;
mRescaledPhyPaths=cell(1,size(cOLPaths,2));
orienScale=2*pi/numOriens;
for j=1:size(cOLPaths,2)
    kappa=CurvatureEstimationFromPath(cOLPaths{j},smoothWidth);
    turningAngle=mod(cOLPaths{j}(3,:),2*pi);
    mRescaledPhyPath=RescaledCoords(cOLPaths{j}(1:2,:),origin,[gridScale;gridScale]);
    mRescaledPhyPaths{j}=mRescaledPhyPath;
    [curvatureMap,scoreMapTem_c]=EstimatePathCurvature(mRescaledPhyPath,kappa,options.imageSize);
    [orienMap,scoreMapTem_o]=EstimatePathCurvature(mRescaledPhyPath,turningAngle,options.imageSize);
    scoreMap_c=scoreMap_c+scoreMapTem_c;
    scoreMap_o=scoreMap_o+scoreMapTem_o;
    curvature_phi=curvature_phi+curvatureMap;
    orien_phi=orien_phi+orienMap;
end

mLinearCoords=find(scoreMap_c>0.5);clear scoreMap_c;
initialVelocity=curvature_phi(mLinearCoords);
initialOriens=orien_phi(mLinearCoords);
[mRows,mCols]=ind2sub(size(curvature_phi),mLinearCoords);
mSeeds=[mCols';mRows'];
extendedVelocityPhysical_Curvature = ComputeExtendedVelocity(mSeeds,initialVelocity,options.imageSize,mask);
extendedVelocity_Orien = ComputeExtendedVelocity(mSeeds,initialOriens,options.imageSize,mask);
extendedCurvatureVelocity=zeros(numRows,numCols,numOriens);
for i=1:numRows
    for j=1:numCols
        if mask(i,j)<0.5
            continue;
        end
        turningAngle=mod(extendedVelocity_Orien(i,j),2*pi);
        turningAngleIndex=mod(round(turningAngle./orienScale)+1,numOriens);
        if turningAngleIndex==0
            turningAngleIndex=numOriens;
        end
        curvatureValue=extendedVelocityPhysical_Curvature(i,j);
        extendedCurvatureVelocity(i,j,:)=curvatureValue;
        extendedCurvatureVelocity(i,j,turningAngleIndex)=-curvatureValue;
        refAngle=(turningAngleIndex-1)*orienScale;
        for m=1:round(numOriens)
            currentAngleindex=mod(turningAngleIndex+m,numOriens); 
            if currentAngleindex==0
                currentAngleindex=numOriens;
            end
            currentAngle=(currentAngleindex-1)*orienScale;
            scalarProduct=cos(refAngle)*cos(currentAngle)+sin(refAngle)*sin(currentAngle);
            extendedCurvatureVelocity(i,j,currentAngleindex)=-curvatureValue*sign(scalarProduct);
        end
    end
end

phyCurvatureMaps.extendedVelocity_Orien=extendedVelocity_Orien;
phyCurvatureMaps.extendedVelocityPhysical_Curvature=extendedVelocityPhysical_Curvature;
phyCurvatureMaps.curvature_phi=curvature_phi;
end

function [curvatureMap,scoreMap]=EstimatePathCurvature(rescaledPath,c,imageSize)
numRows=imageSize(1);
numCols=imageSize(2);
curvatureMap=zeros(numRows,numCols);
scoreMap=zeros(numRows,numCols);

for j=1:size(rescaledPath,2)
    mPx=round(rescaledPath(1,j));
    mPy=round(rescaledPath(2,j));
    curvatureMap(mPy,mPx)=curvatureMap(mPy,mPx)+c(j);
    scoreMap(mPy,mPx)=scoreMap(mPy,mPx)+1;
end

curvatureMap=curvatureMap./(scoreMap+eps);

end

function extendedVelocity = ComputeExtendedVelocity(mSeeds,initialVelocity,imageSize,mask)

numRows=imageSize(1);
numCols=imageSize(2);

input.model='Isotropic2';
input.order=1.0;
input.seedRadius=0;
input.eps=0.1;
input.dims = [numCols;numRows];
input.arrayOrdering='YXZ_ColumnMajor';
input.verbosity=0;
input.origin=[0.0;0.0];  % Physical origin
input.gridScale=1.0; % Physical gridScale
input.exportValues=1.0;
input.exportExtendedVelocity=1.0;
input.initialVelocity=initialVelocity;
input.seeds=mSeeds-1.0;
input.cost=1.0;
input.mask=mask;
output=HamiltonVelocityExtension(input);
extendedVelocity=output.extendedVelocity;
end
