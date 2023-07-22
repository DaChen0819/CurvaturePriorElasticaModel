function c1=CurvatureEstimationFromPath(cOLPath,width)

assert(size(cOLPath,1)==3);
if nargin==1
    width=16;
end

p=cOLPath(1,1:end);
q=cOLPath(2,1:end);
o=cOLPath(3,1:end);
gp=gradient(p,1);
gq=gradient(q,1);
go=gradient(o,1);
c=go./sqrt(gp.^2+gq.^2);
c1=smoothdata(c,'rlowess',width);

end