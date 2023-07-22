function score=ComputeJaccardScore(A,B)

assert(nargin==2);
A=double(A);
A=A./max(A(:));

B=double(B);
B=B./max(B(:));

A=double(A>0.5);
B=double(B>0.5);

cup=double((A+B)>0.5);
cap=double((A.*B)>0.5);
score=sum(cap(:))/sum(cup(:));

end