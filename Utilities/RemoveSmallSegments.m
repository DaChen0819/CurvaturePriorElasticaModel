function newBW=RemoveSmallSegments(binaryImg,thresh)

binaryImg=binaryImg>0.5;
CC = bwconncomp(binaryImg);
newBW=double(binaryImg);
NumObjects=CC.NumObjects;
PixelIdxList=CC.PixelIdxList;

for i=1:NumObjects
    PixelIdx=PixelIdxList{i};
    if size(PixelIdx,1)<=thresh
        newBW(PixelIdx)=0.;
    end
end

end