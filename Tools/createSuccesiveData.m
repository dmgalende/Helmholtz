clear all
close all

%---------------------------------
initialData = 'data_1_1.mat';
meshFiles = {
             'mesh_0.75_1.dcm'
             'mesh_0.75_2.dcm'
%              'mesh_0.25_1.dcm'
%              'mesh_0.25_2.dcm'
%              'mesh_3_1.dcm'
             };
% meshFiles = {};
% for i = 1:3
%     meshFiles = {meshFiles{:} ['Scat_1_' num2str(i) '.dcm']};
% end
% meshFiles = meshFiles';

bottomFiles = {
               'ellipticShoal_aux.m'
               };
smoothBottom = false;
           
outputName = {
             'data_0.75_1.mat'
             'data_0.75_2.mat'
%              'data_0.25_1.mat'
%              'data_0.25_2.mat'
%              'data_3_1.mat'
              };
% outputName = {};
% for i = 1:3
%     outputName = {outputName{:} ['Scat_1_' num2str(i) '.dcm']};
% end
% outputName = outputName';
%---------------------------------

currentDirectory = pwd;
cd ..
path1 = [pwd '/Callbacks'];
path2 = [pwd '/FEM'];
path3 = [path2 '/matlabShapeFunctions'];
path4 = [path3 '/orthopoly'];
addpath(path1,path2,path3,path4)
cd(currentDirectory)

nMesh = length(meshFiles);
for imesh = 1:nMesh
    
    % Initial data
    load(initialData)
    data.plotBottom = [];
    data.plotBoundary = [];
    data.plotMesh = [];
    nOfBoundaries = numel(data.mesh.boundaryNames);
    data.plotPMLarea = cell(2,nOfBoundaries);
    data.plotPMLnodes = cell(2,nOfBoundaries);
    data.PML = cell(6,nOfBoundaries);
    data.PML(5,:) = {'off'};
    data.plotSubmesh = [];
    meshfile = meshFiles{imesh};
    bottomfile = bottomFiles{1}; %bottomFiles{imesh};
    
    % Mesh & Bottom
    data = createMeshDataAUX(data,meshfile);
    data = createBottomDataAUX(data,bottomfile);
    if strcmp(data.computation,'NEFEM')
        nameNodal = 'X';
        nameFac2d = 'faceNodes';
        nameCon = 'T';
        nameType = 'type';
        nameNodes = 'nOfNodes';
        nameElem = 'elemInfo';
    else
        nameNodal = data.mesh.fieldNames{data.mesh.indexElemPosCon(2)};
        nameFac2d = data.mesh.elemFieldNames{data.mesh.indexTypNodFac1Fac2(4)};
        nameCon = data.mesh.fieldNames{data.mesh.indexElemPosCon(3)};
        nameType = data.mesh.elemFieldNames{data.mesh.indexTypNodFac1Fac2(1)};
        nameNodes = data.mesh.elemFieldNames{data.mesh.indexTypNodFac1Fac2(2)};
        nameElem = data.mesh.fieldNames{data.mesh.indexElemPosCon(1)};
    end
    X = data.mesh.(nameNodal);
    T = data.mesh.(nameCon);
    nOfNodesPerElement = size(T,2);
    data.plotCondition = size(X,1) < 1.5e6 && size(T,1) < 4e5;

    % Reference element
    data.mesh.referenceElement = createReferenceElement(data.mesh.(nameElem).(nameType),...
        data.mesh.(nameElem).(nameNodes),[]);

    % PML
    for iboundary = 1:nOfBoundaries
        if data.BC.values(iboundary) == 4
            data.currentBoundary.value = iboundary;
            data.PML{5,data.currentBoundary.value} = 'on';
            nameBoundary = data.mesh.fieldNames{data.mesh.boundaryIndex(data.currentBoundary.value)};
            if any(strcmp(data.computation,{'FEM' 'CDG' 'DG'}))
                    [data.PML{2,data.currentBoundary.value} data.PML{3,data.currentBoundary.value}] = ...
                    getPMLboundaries(data.mesh.(nameNodal),data.mesh.(nameBoundary),...
                                     data.mesh.referenceElement,...
                                     []);
                elseif strcmp(data.computation,'NEFEM')
                    [data.PML{2,data.currentBoundary.value} data.PML{3,data.currentBoundary.value}] = ...
                    getPMLboundaries(data.mesh.X,data.mesh.(nameBoundary),...
                                     data.mesh.referenceElement,...
                                     [],...
                                     data.mesh.nurbs,...
                                     data.mesh.trimmedInfo.(nameBoundary));
            end

            param = data.BC.parameters{data.currentBoundary.value}{4};
            if param(4) > 1
                pos = param(4) - 1;
                nameSubCon = data.mesh.fieldNames{data.mesh.submeshIndex(pos)};
                PMLnodes = getNodes(data.mesh.(nameSubCon));
                data.PML{6,data.currentBoundary.value} = PMLnodes;
            else
                PMLnodes = 1:size(data.mesh.(nameNodal),1);
                data.PML(6,data.currentBoundary.value) = {[]};
            end
            [sigma,sigmaPos,stretching] = coefPML(data.mesh.(nameNodal),PMLnodes,...
                                   data.PML{2,data.currentBoundary.value},...
                                   param(1),param(2),param(3));
            data.PML{1,data.currentBoundary.value} = zeros(size(data.mesh.(nameNodal),1),2);
            data.PML{1,data.currentBoundary.value}(PMLnodes,:) = sigma;
            data.PML{4,data.currentBoundary.value} = sigmaPos;
            data.PML{7,data.currentBoundary.value} = stretching;
        end
    end
    
    %Smooth bottom in PML
    if smoothBottom
        data = createBottomDataAUX(data,0);
    end

    % Save data
    save(outputName{imesh},'data')
    disp(['Saved ' outputName{imesh}])
    clear data X T
end


    




