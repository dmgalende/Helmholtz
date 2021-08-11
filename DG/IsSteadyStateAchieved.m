function [iCycle,errSteadyState,boolSteadyState,refSol] = IsSteadyStateAchieved...
    (iCycle,errSteadyState,theReferenceElement,...
    nOfCyclesMax,U,refSol,theMesh,tolSteadyState,infoFaces,sstBoundary)

X = theMesh.X;
T = theMesh.T;
infoFaces_SST = infoFaces.(['exteriorFaces_',sstBoundary]);
boolSteadyState = 0;
nOfCheckFaces = size(infoFaces_SST,1);
nOfFaceNodes = numel(theReferenceElement.NodesCoord1d);
iCycle = iCycle+1;
if iCycle == 1
    errSteady = 1;
    refSol=zeros(nOfFaceNodes,nOfCheckFaces);
    for i = 1:nOfCheckFaces
        infoFace = infoFaces_SST(i,:);
        elem = infoFace(1);
        face = infoFace(2);
        faceNodes = theReferenceElement.faceNodes(face,:);
        refSol(:,i) = U(faceNodes,4,elem);
    end
else
    refSolPrev = refSol;
    refSol=zeros(nOfFaceNodes,nOfCheckFaces);
    err = zeros(nOfCheckFaces,1);
    solNorm = err;
    for i = 1:nOfCheckFaces
        infoFace = infoFaces_SST(i,:);
        elem = infoFace(1);
        face = infoFace(2);
        faceNodes = theReferenceElement.faceNodes(face,:);
        Te = T(elem,:);
        Xe = X(Te,:);
        Xf = Xe(faceNodes,:);
        refSol(:,i) = U(faceNodes,4,elem);
        err(i) = computeL2NormBoundary(theReferenceElement,Xf,1:nOfFaceNodes,refSol(:,i)-refSolPrev(:,i));
        solNorm(i) = computeL2NormBoundary(theReferenceElement,Xf,1:nOfFaceNodes,refSol(:,i));
    end
    errSteady = sqrt(sum(err.^2)/sum(solNorm.^2));
end
errSteadyState(iCycle)=errSteady;
disp(sprintf('-----------Cycle=%d --- Error Steady=%e -',iCycle,errSteady));
if errSteady<tolSteadyState
    disp(sprintf('Steady state reached after %d cycles', iCycle));
    boolSteadyState = 1;
elseif iCycle==nOfCyclesMax
    disp(sprintf('Steady state not reached after %d cycles', iCycle))
end



