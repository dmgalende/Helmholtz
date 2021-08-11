function faceElmentInfo = computeElementFaceInfo(Tb,N,T,elems_flag)

faces = [1 2 ; 2 3 ; 3 1];
nOfElements = size(Tb,1);
faceElmentInfo = zeros(nOfElements,2);

for i = 1:nOfElements
    vnode1 = Tb(i,1);
    vnode2 = Tb(i,end);
    elem1 = N(vnode1,logical(N(vnode1,:)));
    elem2 = N(vnode2,logical(N(vnode2,:)));
    elem_aux = intersect(elem1,elem2);
    elem = [];
    for j = 1:length(elem_aux)
        pos = elem_aux(j) == elems_flag;
        elems_flag_aux = elems_flag(pos);
        if isscalar(elems_flag_aux), elem = elem_aux(j); break, end
    end
    if isempty(elem), error('error in the ELEMENT computeElementFaceInfo.m'); end
    
    vnode = [vnode1,vnode2];
    if all(vnode == T(elem,faces(1,:)))
        face = 1;
    elseif all(vnode == T(elem,faces(2,:)))
        face = 2;
    elseif all(vnode == T(elem,faces(3,:)))
        face = 3;
    else
        error('error in the FACE computeElementFaceInfo.m');
    end
    
    faceElmentInfo(i,:) = [elem,face];
end
    