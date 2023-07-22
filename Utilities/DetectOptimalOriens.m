function orienIdx=DetectOptimalOriens(mPoints,cost)

orienIdx=[];

for j=1:size(mPoints,2)
    listOS=squeeze(cost(mPoints(2,j),mPoints(1,j),:));
    [~,I]=max(listOS(:));
    orienIdx=cat(2,orienIdx,I(1));
end

end