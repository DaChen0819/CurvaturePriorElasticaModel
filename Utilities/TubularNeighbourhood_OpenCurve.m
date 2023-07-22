function tubularNeighRegion=TubularNeighbourhood_OpenCurve(cOLPaths,options)

gridScale=options.gridScale;
origin=options.origin;
numRows=options.imageSize(1);
numCols=options.imageSize(2);
maxRadius=options.maximalRadius;
seedFlags=[];
mSeeds=[];
mFirstExtendedPoints=[];
mLastExtendedPoints=[];
mFirstPoints=[];
mLastPoints=[];
mDigitalPaths=cell(1,size(cOLPaths,2));
for i=1:size(cOLPaths,2)
    cOLPath=cOLPaths{i};
    if isempty(cOLPath)
        continue;
    end
    cPath=cOLPath(1:2,:);
    mPath=RescaledCoords(cPath,origin,[gridScale;gridScale]);
    [a,mDigitalPath]=ConvertToDigitalPath(mPath,numRows,numCols);
    mDigitalPaths{i}=mDigitalPath;
    [mRows,mCols]=find(a>0);
    seedFlags=cat(2,seedFlags,i.*ones(1,size(mRows,1)));
    mSeeds=cat(2,mSeeds,[mCols';mRows']);
end
voronoiMap=ConstructVoronoiDiagram(options,maxRadius,mSeeds,seedFlags);
tubularNeighRegion=zeros(numRows,numCols);
for i=1:size(cOLPaths,2)
    cOLPath=cOLPaths{i};
    if isempty(cOLPath)
        continue;
    end
    firstOrien=cOLPath(3,end)+pi;
    lastOrien=cOLPath(3,1);
    firstTangent=[cos(firstOrien);sin(firstOrien)];
    lastTangent=[cos(lastOrien);sin(lastOrien)];
    mDigitalPath=mDigitalPaths{i};
    if isempty(mDigitalPath)
        continue;
    end
    mFirstPoint=mDigitalPath(:,end);
    mLastPoint=mDigitalPath(:,1);
    for j=1:0.3:100
        mFirstExtendedPoint=round(j.*firstTangent+mFirstPoint);
        if mFirstExtendedPoint(1)<=0 || mFirstExtendedPoint(2)<=0 || mFirstExtendedPoint(2)>numRows || mFirstExtendedPoint(1)>numCols
            mFirstExtendedPoint=mFirstPoint;
            break; 
        elseif mFirstExtendedPoint(1)~=mFirstPoint(1) || mFirstExtendedPoint(2)~=mFirstPoint(2)
            break;
        end
    end
    for j=1:0.3:100
        mLastExtendedPoint=round(j.*lastTangent+mLastPoint);
        if mLastExtendedPoint(1)<=0 || mLastExtendedPoint(2)<=0 || mLastExtendedPoint(2)>numRows || mLastExtendedPoint(1)>numCols
            mLastExtendedPoint=mLastPoint;
            break; 
        elseif mLastExtendedPoint(1)~=mLastPoint(1) || mLastExtendedPoint(2)~=mLastPoint(2)
            break;   
        end
    end
    mFirstPoints=cat(2,mFirstPoints,mFirstPoint);
    mLastPoints=cat(2,mLastPoints,mLastPoint);
    mFirstExtendedPoints=cat(2,mFirstExtendedPoints,mFirstExtendedPoint);
    mLastExtendedPoints=cat(2,mLastExtendedPoints,mLastExtendedPoint);
    options.mask=double(round(voronoiMap)==i);
    distance1=ComputeGeodesicDistances(mDigitalPath,options);
    distance2=ComputeGeodesicDistances(mFirstExtendedPoint,options);
    distance3=ComputeGeodesicDistances(mLastExtendedPoint,options);
    tubularNeighRegion=tubularNeighRegion+i.*double(and(distance1<distance2,distance1<distance3));
end

end

function distance=ComputeGeodesicDistances(mSeeds,options)

numRows=options.imageSize(1);      
numCols=options.imageSize(2);

input.dims=[numCols;numRows];
input.gridScale=1.0;
input.origin=[0.0;0.0];
input.exportValues=1; 

input.factoringMethod='None';
input.verbosity=0;
input.arrayOrdering='YXZ_ColumnMajor';
input.order=1;
input.exportGeodesicFlow=0.0;
input.exportEuclideanLengths=0.0;
input.model='Isotropic2';
input.cost=1.0;
input.seeds=mSeeds-1.0;
input.stopAtDistance=options.maximalRadius;
input.mask=options.mask;
output_GeoDistance=HamiltonGeodesicDistanceComputation(input);

distance=output_GeoDistance.values;
end

function voronoiMap=ConstructVoronoiDiagram(options,maxRadius,mSeeds,seedFlags)

numRows=options.imageSize(1);
numCols=options.imageSize(2);

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
input.seeds=mSeeds-1.0;
input.seedFlags=seedFlags';
input.exportEuclideanLengths=0.0;

output=HamiltonGeodesicDistanceComputation(input);
voronoiMap=output.voronoiFlags;


end


