function [XX,XY,YX,YY] = faceContribution(theReferenceElement,Xe1,Xe2,n,A,infoFace,...
    alpha_f,C11,C12,sigma_f,omega,Nx_face,Ny_face)

nOfElementNodes = size(theReferenceElement.NodesCoord,1);
nOfFaceNodes = numel(theReferenceElement.faceNodes1d);
iFace =infoFace(2);
iFace2 =infoFace(4);
IPw = theReferenceElement.IPweights1d;
N1 = theReferenceElement.N1d;
N2 = flipud(N1);
ngauss = numel(IPw);

% shape functions
N1xi = Nx_face(:,:,iFace);
N2xi = Nx_face(:,:,iFace2);
N1eta = Ny_face(:,:,iFace);
N2eta = Ny_face(:,:,iFace2);

% x and y coordinates of the element nodes
xe1 = Xe1(:,1); ye1 = Xe1(:,2);
xe2 = Xe2(:,1); ye2 = Xe2(:,2);
%inizialization
XX = zeros(nOfElementNodes);
XY = zeros(nOfElementNodes);
YX = zeros(nOfElementNodes);
YY = zeros(nOfElementNodes);
Con3_11 = zeros(nOfElementNodes,nOfFaceNodes);
Con3_12 = Con3_11; Con3_21 = Con3_11; Con3_22 = Con3_11;
Con7_11 = zeros(nOfFaceNodes,nOfFaceNodes); Con7_12 = Con7_11;
for g = 1:ngauss
    %Shape functions at the current integration point
    N1_g = N1(g,:);
    N2_g = N2(g,:);
    %shape functions derivatives for the actual element
    N1xi_g = N1xi(g,:);
    N1eta_g = N1eta(g,:);
    %shape functions derivatives for the neighbouring element
    N2xi_g = N2xi(ngauss - g+1,:);
    N2eta_g = N2eta(ngauss - g+1,:);
    %bottom parameters at current integration point
    alpha_g = N1_g*alpha_f;
    %PML parameters
    sigma_g = N1_g*sigma_f;
    param = 1 + sqrt(-1)*sigma_g/omega;
    %Jacobian 1
    J = [N1xi_g*xe1	  N1xi_g*ye1
        N1eta_g*xe1  N1eta_g*ye1];
    invJ = inv(J);
    %Jacobian 2
    J2 = [N2xi_g*xe2	  N2xi_g*ye2
        N2eta_g*xe2  N2eta_g*ye2];
    invJ2 = inv(J2);
    %shape functions derivatives for the neighbouring element
    N1x_g = invJ(1,1)*N1xi_g + invJ(1,2)*N1eta_g;
    N1y_g = invJ(2,1)*N1xi_g + invJ(2,2)*N1eta_g;
    N2x_g = invJ2(1,1)*N2xi_g + invJ2(1,2)*N2eta_g;
    N2y_g = invJ2(2,1)*N2xi_g + invJ2(2,2)*N2eta_g;
    %Integration weight
    dL = A/2*IPw(g);
    P_gradNdotn = param(2)/param(1)*N1x_g * n(1) + param(1)/param(2)*N1y_g * n(2);
    P_gradNdotn2 = param(2)/param(1)*N2x_g * n(1) + param(1)/param(2)*N2y_g * n(2);

    Con3_11 = Con3_11 + 0.5 * transpose(P_gradNdotn) * N1_g * alpha_g * dL;
    Con3_12 = Con3_12 + 0.5 * transpose(P_gradNdotn) * N2_g * alpha_g * dL;
    Con3_21 = Con3_21 - 0.5 * transpose(P_gradNdotn2) * N1_g * alpha_g * dL; % il meno è dovuto alla normale..
    Con3_22 = Con3_22 - 0.5 * transpose(P_gradNdotn2) * N2_g * alpha_g * dL; % il meno è dovuto alla normale..

    Con7_11 = Con7_11 +  N1_g' * N1_g * C11 * dL;
    Con7_12 = Con7_12 +  N1_g' * N2_g * C11 * dL;

end
Con4_11 = transpose(Con3_11);
Con4_12 = -transpose(Con3_21);
Con4_21 = -transpose(Con3_12);
Con4_22 = transpose(Con3_22);
Con5_11 = C12 * Con3_11;
Con5_12 = C12 * Con3_12;
Con5_21 = -C12 * Con3_21;
Con5_22 = -C12 * Con3_22;
Con6_11 = C12 * Con4_11;
Con6_12 = C12 * Con4_12;
Con6_21 = -C12 * Con4_21;
Con6_22 = -C12 * Con4_22;
Con7_21 = Con7_12;
Con7_22 = Con7_11;

% assembly
nodes = theReferenceElement.faceNodes(iFace,:);
nodes2 = theReferenceElement.faceNodes(iFace2,:);

XX(:,nodes) = XX(:,nodes) - (+ Con3_11) -(+ Con5_11);
XX(nodes,:) = XX(nodes,:) -(+ Con4_11) -(+ Con6_11);
XX(nodes,nodes) = XX(nodes,nodes) +(+ Con7_11);

XY(:,nodes2) = XY(:,nodes2) -(- Con3_12) -(- Con5_12);
XY(nodes,:) = XY(nodes,:) -(+ Con4_12) -(- Con6_12);
XY(nodes,nodes2) = XY(nodes,nodes2) +(- Con7_12);

YX(:,nodes) = YX(:,nodes) -(- Con3_21) -(- Con5_21);
YX(nodes2,:) = YX(nodes2,:) -(+ Con4_21) -(- Con6_21);
YX(nodes2,nodes) = YX(nodes2,nodes) +(- Con7_21);

YY(:,nodes2) = YY(:,nodes2) - (+ Con3_22) -(+ Con5_22);
YY(nodes2,:) = YY(nodes2,:) -(+ Con4_22) -(+ Con6_22);
YY(nodes2,nodes2) = YY(nodes2,nodes2) +(+ Con7_22);




