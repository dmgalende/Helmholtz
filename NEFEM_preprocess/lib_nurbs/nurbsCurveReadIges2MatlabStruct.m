function nurbs = nurbsCurveReadIges2MatlabStruct(ParameterData)
%
% nurbs = nurbsCurveReadIges2MatlabStruct(ParameterData)
%

nPos = length(ParameterData);
iNurbs = 0;
for iPos = 1:nPos
    if strcmp(ParameterData{iPos}.name,'B-NURBS CRV')
        
        readNurbs = ParameterData{iPos}.nurbs;
        
        nurbsAux.U = readNurbs.knots;
        nurbsAux.Pw = readNurbs.coefs';
        nurbsAux.iniParam = nurbsAux.U(1);
        nurbsAux.endParam = nurbsAux.U(end);
        
        iNurbs = iNurbs + 1;
        nurbs(iNurbs) = nurbsAux;
    end
end
