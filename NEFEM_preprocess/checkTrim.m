function checkTrim(file,tol)

addpath('lib_nurbs')

%------------------------ Detect boundaries
nefemFile = load(file);
nurbs = nefemFile.nurbs;
X = nefemFile.X;
problemNode = false(size(X,1),1);
meshFields = fieldnames(nefemFile);
k = [];
for i = 1:length(meshFields)
    name = meshFields{i};
    if length(name) > 3 && strcmpi(name(1:3),'Tb_'), k = [k i]; end
end

%------------------------ Check trim
for i = k
    T = nefemFile.(meshFields{i});
    disp(['-----> Connectivity ' meshFields{i}])
    trimInfo = nefemFile.trimmedInfo.(meshFields{i});
    for elem = 1:size(T,1)
        nodes = T(elem,1:2); %first face is the curved one
        Xe = X(nodes,:);
        u1 = trimInfo(elem).trim(1);
        u2 = trimInfo(elem).trim(2);
        iNurb = trimInfo(elem).idNurbs;
        aNurb = nurbs(iNurb);
        pt1 = nurbsCurvePoint(aNurb,u1);
        pt2 = nurbsCurvePoint(aNurb,u2);
        dis1 = norm(pt1(1:2) - Xe(1,:))/norm(Xe(1,:));
        dis2 = norm(pt2(1:2) - Xe(2,:))/norm(Xe(2,:));
        if dis1 > tol && ~problemNode(nodes(1))
            disp(['Problem in node ' num2str(nodes(1)) ' elem ' num2str(elem) ' nurbs ' num2str(iNurb) ' param ' num2str(u1)])
            problemNode(nodes(1)) = true;
        end
        if dis2 > tol && ~problemNode(nodes(2))
            disp(['Problem in node ' num2str(nodes(2)) ' elem ' num2str(elem) ' nurbs ' num2str(iNurb) ' param ' num2str(u2)])
            problemNode(nodes(2)) = true;
        end
    end
end


        
        
        
