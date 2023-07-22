function geoLengths=ComputeSpatialGeodesicDistance(options,mSeeds,geoModel,geoMetric,gridParas)

assert(nargin==2 || nargin==4 || nargin ==5);
if isfield(options,'imageSize')
    numRows=options.imageSize(1);
    numCols=options.imageSize(2);
else
    error('imageSize must be specified.');
end

if nargin==5
    gridScale=gridParas.gridScale;
    origin=gridParas.origin;
else
    gridScale=1.0;
    origin=[0.0;0.0];
end

if nargin==2
    geoModel='Isotropic2';
    geoMetric=ones(numRows,numCols);
end

if isfield(options,'exportGeodesicFlow')
    exportGeodesicFlow=options.exportGeodesicFlow;
else
    exportGeodesicFlow=0.0;
end

if isfield(options,'exportEuclideanLengths')
    exportEuclideanLengths=options.exportEuclideanLengths;
else
    exportEuclideanLengths=0.0;
end

if exportEuclideanLengths>0.5
    if isfield(options,'euclideanScale')
        euclideanScale=options.euclideanScale;
    else
        euclideanScale=ones(2,numRows,numCols);
    end
    input_GeoDistance.euclideanScale=euclideanScale;
end


input_GeoDistance.dims=[numCols;numRows];
input_GeoDistance.gridScale=1.0;
input_GeoDistance.origin=[0.0;0.0];
input_GeoDistance.exportValues=1; % distance table, of size [n,n] (minimum of the previous one over directions)
input_GeoDistance.exportActiveNeighs=0;

input_GeoDistance.factoringMethod='None';
input_GeoDistance.showProgress=0;
input_GeoDistance.verbosity=0;
input_GeoDistance.arrayOrdering='YXZ_ColumnMajor';
input_GeoDistance.cosAngleMin=0.5;
input_GeoDistance.refineStencilAtWallBoundary=0;
input_GeoDistance.order=1;
input_GeoDistance.spreadSeeds=-1;
input_GeoDistance.exportGeodesicFlow=exportGeodesicFlow;
input_GeoDistance.exportEuclideanLengths=exportEuclideanLengths;
input_GeoDistance.model=geoModel;

if strcmp(geoModel,'Isotropic2')
    input_GeoDistance.cost=geoMetric;
else
    input_GeoDistance.metric=geoMetric;
end

input_GeoDistance.seeds= RescaledCoords(mSeeds,origin,[gridScale;gridScale],false);

if isfield(options,'maxActiveRegionWidth')
    input_GeoDistance.stopAtDistance=options.maxActiveRegionWidth;
end
if isfield(options,'mask')
    input_GeoDistance.mask=options.mask;
end

output_GeoDistance=HamiltonGeodesicDistanceComputation(input_GeoDistance);
geoLengths=output_GeoDistance.values;

end


