function [X] = one2twoDmapping(iFace, theReferenceElement)
% function to map the 1D gauss points onto
% the iFace face of the reference triangle
% Reference face is [-1,1]
% Reference triangle is [-1,-1; 1,-1; -1,1]
% counter-clockwise mapping
if iFace == 1
    x = theReferenceElement.IPcoordinates1d;
    y = -1 * ones(size(x));
elseif iFace ==2
    x = -theReferenceElement.IPcoordinates1d;
    y = -x;
elseif iFace ==3
    y = -theReferenceElement.IPcoordinates1d;
    x = -1 * ones(size(y));
end
X =[x,y];
