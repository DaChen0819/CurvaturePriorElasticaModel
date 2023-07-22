function mPoints=FindClickedPoints(rawImg,numSeeds)
 if numSeeds<=0
     mPoints=[];
     return;
 end
figure(1);
imagesc(rawImg); axis off;axis equal;
title('Click Points');
if size(rawImg,3)~=3
    colormap(gray);
end
hold on;
mPoints=zeros(2,numSeeds);
for i=1:numSeeds
    [mPoints(1,i),mPoints(2,i)] = ginput( 1 );
    plot(mPoints(1,i),mPoints(2,i), 'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',10);
    pause(0.001)
end
hold off
mPoints=round(mPoints);

pause(1);
close all;

end
