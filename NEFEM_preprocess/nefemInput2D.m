%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script creates the required mesh for a NEFEM computation using the
% Berkhoff GUI. The script assumes that the only curved elements are those
% that have a UNIQUE face on the boundary. The interior elements are 
% considered as non-curved ones.
%
% Required inputs:
%       1 - .igs file with the exact geometry.
%       3 - .dcm EZ4U file with the mesh of an approximate geometry. This
%           file has to contain all the boundaries as attributes sorted by
%           elements (directly choosen in EZ4U).
%
% Output: .mat file containing...
%       1 - nodal position matrix (X).
%       2 - 2D interior connectivity matrix (T).
%       3 - 2D boundary connectivity matrices (Tb_(boundaryName)). The first
%           face always belongs to the boundary.
%       4 - a struct (nurbs) with the exact geometry.
%       5 - a struct (trimmedInfo) with the following fields for each
%           boundary and each element:    
%           5.1 - idNurbs: corresponding nurb for the element.
%           5.2 - trim: initial and ending parameter for the idNurb in the 
%                       element.
%
% IMPORTANT REQUIRED: there have to be a node on the initial and ending
%                     point of all exact geometry nurbs. The selected 
%                     boundaries from EZ4U have to begin and to end on the
%                     initial and ending point of any exact geometry nurb.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
home

addpath(genpath(pwd))

global nurbs referenceElement referenceElementTri3 tol tolParam projectNodes adaptedNodesToNurbsBoundary

%------------------------ User data
exactGeoIgesFileName = 'ex.igs';
meshFileName         = 'ex_p1.dcm';
outputFileName       = 'ex_p1_nefemMesh_2.mat';

%------------------------ Projection parameter and tolerances in point inversion
tol = 1e-4;
tolParam = 1e-4;
projectNodes = false;
adaptedNodesToNurbsBoundary = false;

%------------------------ Read geometries and EZ4U mesh
nurbs = nurbsCurveReadIges2MatlabStruct(iges2matlab(exactGeoIgesFileName));
if strcmpi(meshFileName(end-3:end),'.mat')
    mesh = load(meshFileName);
else
    meshFileName2Load = GenerateMatFileFromEZ4U_NEFEM(meshFileName);
    mesh = load(meshFileName2Load{:});
    delete(meshFileName2Load{:})
end
referenceElementTri3 = createReferenceElement(1,3);
referenceElement = createReferenceElement(mesh.elemInfo.type,mesh.elemInfo.nOfNodes);
meshFields = fieldnames(mesh);
k = [];
for i = 1:length(meshFields)
    name = meshFields{i};
    if length(name) > 3 && strcmpi(name(1:3),'Tb_'), k = [k i]; end
end

%------------------------ Plot nurbs
nOfNurbs = numel(nurbs);
hold on
for inurb = 1:nOfNurbs
    aNurb = nurbs(inurb);
    h = nurbsCurvePlot(aNurb,aNurb.iniParam,aNurb.endParam,1000);
    set(h,'color','r','linewidth',2)
    pt = nurbsCurvePoint(aNurb,aNurb.iniParam);
    hold on, text(pt(1),pt(2),num2str(inurb))
end
hold off
axis equal
pause(0.1)

%------------------------ Curved and non-curved elements
elem2delete = [];
for i = k
    curvedElem = mesh.(meshFields{i})(:,1);
    curvedFace = mesh.(meshFields{i})(:,2);
    mesh = rmfield(mesh,meshFields{i});
    meshFields{i} = ['Tb_' meshFields{i}(4:end)];
    mesh.(meshFields{i}) = [mesh.T(curvedElem,:) curvedFace];
    elem2delete = [elem2delete ; curvedElem];
end
mesh.T(elem2delete,:) = [];
 
%------------------------ Projecting the curved element vertices to exact geometry and trim
[mesh,trimmedInfo] = projectAndTrimBoundaryNurbs(mesh,meshFields,k);

%------------------------ Plot nurbs tangent vector
hold on
for ivec = 1:nOfNurbs
    aNurb = nurbs(ivec);
    pt = nurbsCurvePoint(aNurb,aNurb.iniParam);
    tnurb = nurbsCurveDerivPoint(aNurb,aNurb.iniParam);
    tnurb = tnurb / norm(tnurb);
    quiver(pt(1),pt(2),tnurb(1),tnurb(2),'color','g','lineWidth',2)
end
hold off

%------------------------ Save
save(outputFileName,'nurbs','trimmedInfo')
save(outputFileName,'-struct','mesh','-append')



