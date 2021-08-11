nPuntos = size(X,1);
[nElementos,nNodosPorElemento] = size(T);

fid = fopen('BarErrorP3Static.vtk','w','d');

fprintf(fid,'# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS %i float\n',nPuntos);

fprintf(fid,'%e %e %e\n',[X zeros(nPuntos,1)]');

fprintf(fid,'\nPOLYGONS %i %i\n',nElementos,nElementos*(nNodosPorElemento+1));

fprintf(fid,'%i %i %i %i\n',[3*ones(nElementos,1) T-1]');

fprintf(fid,'\nPOINT_DATA %i\nFIELD FieldData 1\nCampo 1 %i float\n',nPuntos,nPuntos);

fprintf(fid,'%f\n',campo);

fclose(fid);