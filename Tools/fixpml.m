clear all
clc
dirList =  dir('Barcelona_h10_P*');
meshList = {dirList.name};

% translation
x_trans = 0;
y_trans = 0;

for imesh = 1:numel(meshList)
   meshName = meshList{imesh};
   load(meshName)
   disp(['Loaded ' meshName])
   lim_x_inf = -6560;
   lim_x_sup = -6440;
   lim_y_inf = -865.95;
   lim_y_sup = -865.85;
      X(:,1) = X(:,1)-x_trans;
   X(:,2) = X(:,2)-y_trans;
   x = X(:,1);
   y = X(:,2);

   nodes_1 = Tb_PML(:,1);
   nodes_2 = Tb_PML(:,end);

   aux = x(nodes_1)>lim_x_inf & x(nodes_2)>lim_x_inf & x(nodes_1)<lim_x_sup & x(nodes_2)<lim_x_sup...
       & y(nodes_1)>lim_y_inf & y(nodes_2)>lim_y_inf & y(nodes_1)<lim_y_sup & y(nodes_2)<lim_y_sup;

%     plotMesh(X,T);

   nodes = unique(Tb_PML(aux,:));
%     hold on
%     plot(X(nodes,1),X(nodes,2),'ro')

   Tb_PML_2  = Tb_PML(aux,:);
   Tb_PML(aux,:) = [];

   save(meshName,'X','T','Tb*','elemInfo','elementFaceInfo')
   clear X T Tb* elemInfo elementFaceInfo
end 