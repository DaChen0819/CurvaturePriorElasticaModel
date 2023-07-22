function NewCoordinates = RescaledCoords(oldCoordinates,origin,spacing,direct)

    assert(nargin==3 || nargin==4);
    assert(size(origin,2)==1);
    assert(size(oldCoordinates,1)==size(origin,1));
    assert(all(size(origin)==size(spacing)));

    if nargin==3
        direct=true; % form Physical space to MATLAB index system.
    end

    numPoints = size(oldCoordinates,2); 
    origins=repmat(origin,[1,numPoints]);
    spacings=repmat(spacing,[1,numPoints]);

    if direct
        NewCoordinates = 1+(oldCoordinates-origins)./spacings;
    else
        NewCoordinates = origins+(oldCoordinates-1).*spacings;
    end
end
