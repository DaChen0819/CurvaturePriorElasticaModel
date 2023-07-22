clear variables; close all; clc;

%parameters setting.
numOriens=72;
alpha=4;
maximalRadius_Rectified=2;
maximalRadius_Tube=5;
gridScale=2*pi/numOriens;
origin=[0;0];

% load and demonstrate images.
load('data/rawImg.mat');
[numRows,numCols,numChannels]=size(rawImg);
figure(1);imagesc(rawImg);axis off;axis image;colormap gray;

% load the tubular score map used for generating the binary segmentation
% result.
load('data/tubularScore.mat');
figure(2);imagesc(tubularScore);axis off;axis image;colorbar;

% set the thresholding value for binarizing processing.
thresValue=0.35;
binarySeg = tubularScore>thresValue; % curvilinear strcuture segmentation.
skeletonMap = bwskel(binarySeg); % perform the skeletonization processing.
branchPtMap=bwmorph(skeletonMap,'branchpoints'); % find all the branch points
newSkeletonMap=RemoveNeighPointsFromEnlargedNeighbourhood(double(skeletonMap),double(branchPtMap));
seperateSkeletonMap=RemoveSmallSegments(newSkeletonMap,5);
figure(3);imagesc(seperateSkeletonMap);axis off;axis image;colorbar;

% next step is to perform a smoothing procedure to these disjoint
% skeletons. 
load('data/normalizedOS.mat');% normalizedOS is the normalized orientation score map whose values are ranged in [0,1];
options.eps=0.1;
options.numOriens=numOriens;
options.euclideanScale=[1.0;1.0;0.0];
options.model='Elastica2';
options.gridScale=gridScale;
options.origin=origin;
options.admPercent=0.99;
options.xi=4;
options.imageSize=[numRows;numCols];
options.smoothWidth=5;
options.maximalRadius=maximalRadius_Rectified;

try
[voronoiMap,taggedTrajectory]=ConstructVoronoiTessellation(seperateSkeletonMap,options.maximalRadius);
[cPaths_Rectified,~]=RectifySegmentEndPoints(voronoiMap,taggedTrajectory,exp(-alpha.*normalizedOS),options);
catch
load('cPaths_Rectified.mat');
end


% exhibit the smoothed candidate centerlines.
figure(4);imagesc(ones(numRows,numCols,3));colormap gray;axis image;axis off; hold on;
for j=1:size(cPaths_Rectified,2)
    mRescaledGeo=RescaledCoords(cPaths_Rectified{j}(1:2,:),origin,[gridScale;gridScale]);
    lengthCurve=size(cPaths_Rectified{j},2);
    angle=cPaths_Rectified{j}(3,round(lengthCurve/2));
    plot(mRescaledGeo(1,:),mRescaledGeo(2,:),...
        'color',[0.8*rand(1),min(1.1*rand(1),1),min(1.5*rand(1),0.9)],'LineWidth',3);
    mPt=zeros(2,1);
    mPt(1)=mRescaledGeo(1,round(lengthCurve/2))+6;
    mPt(2)=mRescaledGeo(2,round(lengthCurve/2))+6;
    leng=24;width=3;hatSize=18;
    h3=func_DrawArrow(mPt,mPt+leng*[cos(angle);sin(angle)],hatSize,'BaseAngle',30,'Width',width);
    set(h3,'Facecolor','r','EdgeColor','r');
end

% exhibit the curvature of these candidate centerlines.
smoothWidth=10;
figure(5);imagesc(ones(numRows,numCols,3)); axis image;axis off; hold on;
for j=1:size(cPaths_Rectified,2)   
    mRescaledGeodesic=RescaledCoords(cPaths_Rectified{j}(1:2,:),origin,[gridScale;gridScale]);
    c=CurvatureEstimationFromPath(cPaths_Rectified{j},smoothWidth);
    if j==5 ||j==6
        patch(cat(2,mRescaledGeodesic(1,3:end-1),nan),cat(2,mRescaledGeodesic(2,3:end-1),nan),[c(3:end-1) 0]',...
        'edgecolor','interp','MarkerFaceColor','flat','LineWidth',3);
    else 
        patch(cat(2,mRescaledGeodesic(1,1:end-1),nan),cat(2,mRescaledGeodesic(2,1:end-1),nan),[c(1:end-1) 0]',...
        'edgecolor','interp','MarkerFaceColor','flat','LineWidth',3);
    end

end
h=colorbar('north'); 
colormap turbo;
set(gca,'XTickLabel','');set(gca,'YTickLabel','');




