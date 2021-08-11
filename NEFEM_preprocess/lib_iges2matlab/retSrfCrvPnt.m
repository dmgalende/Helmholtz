function [P,isSCP,isSup,TRI,UV,srfind]=retSrfCrvPnt(SCP,ParameterData,isSup,ind,n,dim)
% RETSRFCRVPNT is a subfunction in IGES2MATLAB file collection.
% No complete documentation is given.
% 
% SCP - 1, surface
%     - 2, curve
%     - 3, point
%     
% ParameterData - Parameter data from IGES file    
%     
% isSup - 1, if superior then return isSup=1
%       - 0, isSup=0 (always)
%       
% ind - index
% 
% n - number of points for curves, n^2 number of poinst for non trimmed surface
% 
% dim - [2,3] 2, curve in domain, 3, curve in space
%                
% m-file can be downloaded for free at
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=
% 13253&objectType=file
% 
% retSrfCrvPnt version 1.3
% 
% written by Per Bergström 2007-01-08
%

isSCP=0;

if SCP==1           % SURFACE
    
    if ParameterData{ind}.type==128
        isSCP=1;
        srfind=ind;
        if and(isSup,ParameterData{ind}.superior)
            P=zeros(3,0);
            TRI=zeros(0,3);
            UV=zeros(2,0);
        else
            
            isSup=0;
            [V,U]=meshgrid(linspace(ParameterData{ind}.v(1),ParameterData{ind}.v(2),n),linspace(ParameterData{ind}.u(1),ParameterData{ind}.u(2),n));
            
            UV=[reshape(U,1,n^2);reshape(V,1,n^2)];
            
            P=nrbevalIGES(ParameterData{ind}.nurbs,UV);
            
            TRI=zeros(2*(n-1)^2,3);
            
            for i=1:(n-1)
                for j=1:(n-1)
                    TRI((2*(i+(j-1)*(n-1))-1),:)=[(i+n*(j-1)) (i+n*(j-1)+1) (i+n*j)];
                    TRI((2*(i+(j-1)*(n-1))),:)=[(i+n*(j-1)+1) (i+n*j) (i+n*j+1)];
                end
            end
        end
        
    elseif ParameterData{ind}.type==144
        
        isSCP=1;
        isSup=0;

        srfind=ParameterData{ind}.pts;

        nO=30;
        NO=4*nO;

        NI=2*nO;
        
        n2=ParameterData{ind}.n2;

        NN=[NO,NI*ones(1,n2)];
        NN=cumsum(NN);

        UV=zeros(2,NN(n2+1));
        
        if ParameterData{srfind}.type==128

            umin=ParameterData{srfind}.u(1);
            umax=ParameterData{srfind}.u(2);
            vmin=ParameterData{srfind}.v(1);
            vmax=ParameterData{srfind}.v(2);

        end

        if ParameterData{ind}.n1
            UV(:,1:NO)=retCrv(ParameterData,ParameterData{ind}.pto,NO,2);
        else
            UV(1,1:nD)=linspace(umin,umax-(umax-umin)/nO,nO);
            UV(2,1:nO)=vmin*ones(1,nO);
            UV(1,(nO+1):(2*nO))=umax*ones(1,nO);
            UV(2,(nO+1):(2*nO))=linspace(vmin,vmax-(vmax-vmin)/nO,nO);
            UV(1,(2*nO+1):(3*nO))=linspace(umax,umin+(umax-umin)/nO,nO);
            UV(2,(2*nO+1):(3*nO))=vmax*ones(1,nO);
            UV(1,(3*nO+1):(4*nO))=umin*ones(1,nO);
            UV(2,(3*nO+1):(4*nO))=linspace(vmax,vmin+(vmax-vmin)/nO,nO);
        end

        if n2>0

            for j=1:n2
                UV(:,((j-1)*NI+1+NO):(j*NI+NO))=retCrv(ParameterData,ParameterData{ind}.pti(j),NI,2);
            end

            ver=char(version);
            vernum=str2num(ver(1:3));

            if vernum>6.99
                TRI = delaunay(UV(1,1:NN(n2+1)),UV(2,1:NN(n2+1)),{'Qt','Qbb','Qc','Pp'});
            else
                warning off
                umax=max(UV(1,1:NN(n2+1)));
                umin=min(UV(1,1:NN(n2+1)));
                vmax=max(UV(2,1:NN(n2+1)));
                vmin=min(UV(2,1:NN(n2+1)));
                dd=sqrt((umax-umin)^2+(vmax-vmin)^2);
                TRI = delaunay(UV(1,1:NN(n2+1)),UV(2,1:NN(n2+1)),dd*1e-4);
                warning on
            end

            [mTRI,nTRI]=size(TRI);

            isInside=logical(ones(mTRI,1));

            for j=1:mTRI

                soTRIj=sort(TRI(j,:));

                if soTRIj(3)<=NO
                    deTRIj=det([(UV(:,soTRIj(2))-UV(:,soTRIj(1))) (UV(:,soTRIj(3))-UV(:,soTRIj(1)))]);
                    if deTRIj<0
                        isInside(j)=0;
                    end

                elseif and(soTRIj(1)>NO,isSameCrv(soTRIj,NN))
                    deTRIj=det([(UV(:,soTRIj(2))-UV(:,soTRIj(1))) (UV(:,soTRIj(3))-UV(:,soTRIj(1)))]);

                    if deTRIj<0
                        isInside(j)=0;
                    end
                end
            end
        else

            ver=char(version);
            vernum=str2num(ver(1:3));            

            if vernum>6.99
                TRI = delaunay(UV(1,1:NO),UV(2,1:NO),{'Qt','Qbb','Qc','Pp'});
            else
                warning off
                umax=max(UV(1,1:NO));
                umin=min(UV(1,1:NO));
                vmax=max(UV(2,1:NO));
                vmin=min(UV(2,1:NO));
                dd=sqrt((umax-umin)^2+(vmax-vmin)^2);               
                TRI = delaunay(UV(1,1:NO),UV(2,1:NO),dd*1e-4);
                warning on
            end

            [mTRI,nTRI]=size(TRI);

            isInside=logical(ones(mTRI,1));

            for j=1:mTRI
                soTRIj=sort(TRI(j,:));
                deTRIj=det([(UV(:,soTRIj(2))-UV(:,soTRIj(1))) (UV(:,soTRIj(3))-UV(:,soTRIj(1)))]);

                if deTRIj<0
                    isInside(j)=0;
                end
            end

        end
        
        if nargin==6
            [UV,TRI]=splitTRI(UV,TRI(isInside,:),sum(isInside));
        else
            TRI=TRI(isInside,:);
        end

        if ParameterData{srfind}.type==128           
            
            P=nrbevalIGES(ParameterData{srfind}.nurbs,UV);
            
        end

    else
        P=zeros(3,0);
        isSCP=0;
        isSup=0;
        TRI=zeros(0,3);
        UV=zeros(2,0);
        srfind=0;
    end

elseif SCP==2       % CURVE
    
    isSCP=0;
    
    if isSup
        if ParameterData{ind}.type==110 | ParameterData{ind}.type==126
            
            isSCP=1;
            
            if ParameterData{ind}.superior
                P=zeros(3,0);
                TRI=0;
            else
                P=retCrv(ParameterData,ind,n,dim);
                isSup=0;
                TRI=0;
            end
        else
            P=zeros(3,0);
            isSup=1;
            TRI=0;            
        end
        
    else
        P=retCrv(ParameterData,ind,n,dim);
        
        if not(isempty(P))
            isSCP=1;
        end
        
        isSup=0;
        TRI=0;
    end

    UV=0;
    srfind=0;
    
elseif SCP==3       % POINT
    
    if ParameterData{ind}.type==116
        P=ParameterData{ind}.p;
        isSCP=1;
        isSup=0;
        TRI=0;
    else
        P=zeros(3,1);
        isSCP=0;
        isSup=1;
        TRI=0;
    end
   
    UV=0;
    srfind=0;
    
end


function rCrv=retCrv(ParameterData,ind,n,dim)

if ParameterData{ind}.type==142

    rCrv=retCrv(ParameterData,ParameterData{ind}.bptr,n,dim);

elseif ParameterData{ind}.type==102

    nvecF=(n/(ParameterData{ind}.length))*(ParameterData{ind}.lengthcnt);
    nvecI=floor(nvecF);
    nrest=n-sum(nvecI);
    nvecre=nvecI-nvecF;

    [so,in]=sort(nvecre);
    nvecI(in(1:nrest))=nvecI(in(1:nrest))+1;

    rCrv=zeros(dim,n);

    stind=1;
    for i=1:(ParameterData{ind}.n)
        if nvecI(i)>0
            endind=stind+nvecI(i)-1;
            rCrv(:,stind:endind)=retCrv(ParameterData,ParameterData{ind}.de(i),nvecI(i),dim);
            stind=endind+1;
        end
    end

elseif ParameterData{ind}.type==110

    if dim==2
        rCrv=zeros(dim,n);
        tvec=linspace(0,(1-1/n),n);
        tvec(2:n)=tvec(2:n)+(rand(1,(n-1))-0.5)/(6*n);
    else
        rCrv=zeros(dim,2);
        tvec=linspace(0,1,2);
    end
    
    rCrv(1,:)=ParameterData{ind}.x1+tvec*(ParameterData{ind}.x2-ParameterData{ind}.x1);
    rCrv(2,:)=ParameterData{ind}.y1+tvec*(ParameterData{ind}.y2-ParameterData{ind}.y1);
    
    if dim>2
        rCrv(3,:)=ParameterData{ind}.z1+tvec*(ParameterData{ind}.z2-ParameterData{ind}.z1);
    end

elseif ParameterData{ind}.type==126

    tvec=zeros(1,n);

    tst=ParameterData{ind}.v(1);
    ten=ParameterData{ind}.v(2);
    
    if dim==2
        tvec=linspace(tst,ten-(ten-tst)/n,n);
        tvec(2:n)=tvec(2:n)+(rand(1,(n-1))-0.5)*((ten-tst)/(6*n));
    else
        tvec=linspace(tst,ten,n);
    end

    p=nrbevalIGES(ParameterData{ind}.nurbs,tvec);

    rCrv=p(1:dim,:);
    
else
    
    rCrv=zeros(dim,0);
    
end


function isSC=isSameCrv(soTRIj,NN)

[ma1,ind1]=max(soTRIj(1)<=NN);
[ma2,ind2]=max(soTRIj(2)<=NN);
[ma3,ind3]=max(soTRIj(3)<=NN);

isSC=and(ind1(1)==ind2(1),ind1(1)==ind3(1));

function [P,Pw]=nrbevalIGES(nurbs,UV)

if size(UV,1)==2

    nUV = size(UV,2);

    if sum(UV(1,:)<nurbs.knots{1}(1))>0
        UV(1,UV(1,:)<nurbs.knots{1}(1))=nurbs.knots{1}(1);
    end
    if sum(UV(1,:)>nurbs.knots{1}(end))>0
        UV(1,UV(1,:)>nurbs.knots{1}(end))=nurbs.knots{1}(end);
    end
    if sum(UV(2,:)<nurbs.knots{2}(1))>0
        UV(2,UV(2,:)<nurbs.knots{2}(1))=nurbs.knots{2}(1);
    end
    if sum(UV(2,:)>nurbs.knots{2}(end))>0
        UV(2,UV(2,:)>nurbs.knots{2}(end))=nurbs.knots{2}(end);
    end

    degree = nurbs.order-1;

    val = reshape(nurbs.coefs,4*nurbs.number(1),nurbs.number(2));
    val = bspeval(degree(2),val,nurbs.knots{2},UV(2,:));
    val = reshape(val,[4 nurbs.number(1) nUV]);

    pnts = zeros(4,nUV);
    for i = 1:nUV
        coefs = squeeze(val(:,:,i));
        pw(:,i) = bspeval(degree(1),coefs,nurbs.knots{1},UV(1,i));
    end

    if nargout==1
        P = pw(1:3,:)./repmat(pw(4,:),3,1);
    else
        P=pw(1:3,:);
        Pw=repmat(pw(4,:),3,1);
    end

elseif size(UV,1)==1
    
    if sum(UV(1,:)<nurbs.knots(1))>0
        UV(1,UV(1,:)<nurbs.knots(1))=nurbs.knots(1);
    end
    if sum(UV(1,:)>nurbs.knots(end))>0
        UV(1,UV(1,:)>nurbs.knots(end))=nurbs.knots(end);
    end    

    val = bspeval(nurbs.order-1,nurbs.coefs,nurbs.knots,UV);
    
    if nargout==1
        P = val(1:3,:)./repmat(val(4,:),3,1);
    else
        P=val(1:3,:);
        Pw=repmat(val(4,:),3,1);
    end

end

function [UV,TRI]=splitTRI(UV,TRI,nIP)

lTRI=size(TRI,1);
lUV=size(UV,2);

UV=[UV,zeros(2,nIP)];

TRIcorner=zeros(lTRI+2*nIP,3);
TRItricorner=zeros(lTRI+2*nIP,3);

TRI=sort(TRI,2);
[TRI,sr]=sortrows(TRI,1:3);

for i=1:(lTRI-1)

    bol=logical(0);

    for j=(i+1):lTRI

        if TRI(j,1)==TRI(i,1)

            if TRI(j,2)==TRI(i,2)

                TRIcorner(i,1)=3;
                TRIcorner(j,1)=3;
                
                TRItricorner(i,1)=j;
                TRItricorner(j,1)=i;
                
            elseif TRI(j,2)==TRI(i,3)

                TRIcorner(i,2)=3;
                TRIcorner(j,1)=2;
                
                TRItricorner(i,2)=j;
                TRItricorner(j,1)=i;                

            elseif TRI(j,3)==TRI(i,3)

                TRIcorner(i,2)=2;
                TRIcorner(j,2)=2;
                
                TRItricorner(i,2)=j;
                TRItricorner(j,2)=i;                

            end
            
        elseif and(TRI(j,2)==TRI(i,2),TRI(j,3)==TRI(i,3))
            
            TRIcorner(i,3)=1;
            TRIcorner(j,3)=1;
            
            TRItricorner(i,3)=j;
            TRItricorner(j,3)=i;            
            
        elseif TRI(j,1)==TRI(i,2)

            if TRI(j,2)==TRI(i,3)

                TRIcorner(i,3)=3;
                TRIcorner(j,1)=1;
                
                TRItricorner(i,3)=j;
                TRItricorner(j,1)=i;
                
                bol=not(bol);

            elseif TRI(j,3)==TRI(i,3)

                TRIcorner(i,3)=2;
                TRIcorner(j,2)=1;
                
                TRItricorner(i,3)=j;
                TRItricorner(j,2)=i;
                
                bol=not(bol);

            end

        elseif TRI(j,1)>TRI(i,2)

            bol=not(bol);            

        end

        if bol
            break
        end

    end
end

TRI=[TRI;zeros(2*nIP,3)];

cTRI=zeros(1,lTRI+2*nIP);

for i=1:lTRI
    cTRI(i)=compTRI(UV(:,TRI(i,:)));
end

for i=1:nIP
    
    cusucTRI=cumsum(cTRI);
    cmp=(rand*0.99999999999999)*cusucTRI(end);
    
    TRIindvec=find(cusucTRI>=cmp);
    TRIind=TRIindvec(1);
   
    UV(:,lUV+i)=triIPoints(UV(:,TRI(TRIind,:)));
    TRI(lTRI+2*(i-1)+1,1:2)=TRI(TRIind,2:3);
    TRI(lTRI+2*(i-1)+2,1:2)=TRI(TRIind,[1,3]);
    
    if TRItricorner(TRIind,2)>0
        TRItricorner(TRItricorner(TRIind,2),(4-TRIcorner(TRIind,2)))=lTRI+2*(i-1)+2;
        TRIcorner(TRItricorner(TRIind,2),(4-TRIcorner(TRIind,2)))=3;
    end
    
    if TRItricorner(TRIind,3)>0
        TRItricorner(TRItricorner(TRIind,3),(4-TRIcorner(TRIind,3)))=lTRI+2*(i-1)+1;
        TRIcorner(TRItricorner(TRIind,3),(4-TRIcorner(TRIind,3)))=3;
    end

    TRI(lTRI+2*(i-1)+1,3)=lUV+i;
    TRIcorner(lTRI+2*(i-1)+1,1)=TRIcorner(TRIind,3);
    TRIcorner(lTRI+2*(i-1)+1,2)=1;
    TRIcorner(lTRI+2*(i-1)+1,3)=1;    
    TRItricorner(lTRI+2*(i-1)+1,1)=TRItricorner(TRIind,3);
    TRItricorner(lTRI+2*(i-1)+1,2)=TRIind;
    TRItricorner(lTRI+2*(i-1)+1,3)=lTRI+2*(i-1)+2;
    
    TRI(lTRI+2*(i-1)+2,3)=lUV+i;
    TRIcorner(lTRI+2*(i-1)+2,1)=TRIcorner(TRIind,2);
    TRIcorner(lTRI+2*(i-1)+2,2)=2;
    TRIcorner(lTRI+2*(i-1)+2,3)=1;
    TRItricorner(lTRI+2*(i-1)+2,1)=TRItricorner(TRIind,2);
    TRItricorner(lTRI+2*(i-1)+2,2)=TRIind;
    TRItricorner(lTRI+2*(i-1)+2,3)=lTRI+2*(i-1)+1;

    TRI(TRIind,3)=lUV+i;
    TRIcorner(TRIind,2)=2;
    TRIcorner(TRIind,3)=2;
    TRItricorner(TRIind,2)=lTRI+2*(i-1)+2;
    TRItricorner(TRIind,3)=lTRI+2*(i-1)+1;

    cTRI(lTRI+2*(i-1)+1)=compTRI(UV(:,TRI(lTRI+2*(i-1)+1,:)));
    cTRI(lTRI+2*(i-1)+2)=compTRI(UV(:,TRI(lTRI+2*(i-1)+2,:)));
    cTRI(TRIind)=compTRI(UV(:,TRI(TRIind,:)));
    
    for j=[(lTRI+2*(i-1)+1),(lTRI+2*(i-1)+2),TRIind]
        Tti=TRItricorner(j,1);
        if Tti>0
            
            Ttii=TRI(Tti,TRIcorner(j,1));
            
            cTRI_I=compTRI(UV(:,[TRI(j,1),Ttii,(lUV+i)]));
            cTRI_II=compTRI(UV(:,[TRI(j,2),Ttii,(lUV+i)]));
            
            if max([cTRI_I,cTRI_II])<5*max(cTRI([j,Tti]))

                if Ttii>TRI(j,2)

                    cTRI(j)=cTRI_I;
                    cTRI(Tti)=cTRI_II;

                    indI=[TRI(j,1) Ttii (lUV+i)];
                    indII=[TRI(j,2) Ttii (lUV+i)];

                    TRIcornerI=[TRIcorner(Tti,2),TRIcorner(j,2),1];
                    TRIcornerII=[TRIcorner(Tti,3),TRIcorner(j,3),1];
                    
                    TRItricornerI=[TRItricorner(Tti,2),TRItricorner(j,2),Tti];
                    TRItricornerII=[TRItricorner(Tti,3),TRItricorner(j,3),j];
                    
                    if TRItricorner(j,3)>0
                        TRItricorner(TRItricorner(j,3),(4-TRIcorner(j,3)))=Tti;
                        TRIcorner(TRItricorner(j,3),(4-TRIcorner(j,3)))=2;
                    end
                    if TRItricorner(Tti,2)>0
                        TRItricorner(TRItricorner(Tti,2),(4-TRIcorner(Tti,2)))=j;
                        TRIcorner(TRItricorner(Tti,2),(4-TRIcorner(Tti,2)))=3;
                    end
                    if TRItricorner(Tti,3)>0
                        TRIcorner(TRItricorner(Tti,3),(4-TRIcorner(Tti,3)))=3;
                    end
 
                    TRI(j,:)=indI;
                    TRI(Tti,:)=indII;
                    
                    TRIcorner(j,:)=TRIcornerI;
                    TRIcorner(Tti,:)=TRIcornerII;
                    
                    TRItricorner(j,:)=TRItricornerI;
                    TRItricorner(Tti,:)=TRItricornerII;                    
                    
                elseif Ttii>TRI(j,1)

                    cTRI(j)=cTRI_I;
                    cTRI(Tti)=cTRI_II;

                    indI=[TRI(j,1) Ttii (lUV+i)];
                    indII=[Ttii TRI(j,2) (lUV+i)];
                    
                    TRIcornerI=[TRIcorner(Tti,1),TRIcorner(j,2),2];
                    TRIcornerII=[TRIcorner(Tti,3),1,TRIcorner(j,3)];
                    
                    TRItricornerI=[TRItricorner(Tti,1),TRItricorner(j,2),Tti];
                    TRItricornerII=[TRItricorner(Tti,3),j,TRItricorner(j,3)];
                    
                    if TRItricorner(j,3)>0
                        TRItricorner(TRItricorner(j,3),(4-TRIcorner(j,3)))=Tti;
                    end
                    if TRItricorner(Tti,1)>0
                        TRItricorner(TRItricorner(Tti,1),(4-TRIcorner(Tti,1)))=j;
                    end
                    if TRItricorner(Tti,3)>0
                        TRIcorner(TRItricorner(Tti,3),(4-TRIcorner(Tti,3)))=3;
                    end

                    TRI(j,:)=indI;
                    TRI(Tti,:)=indII;
                    
                    TRIcorner(j,:)=TRIcornerI;
                    TRIcorner(Tti,:)=TRIcornerII;
                    
                    TRItricorner(j,:)=TRItricornerI;
                    TRItricorner(Tti,:)=TRItricornerII;                    
                    
                else

                    cTRI(j)=cTRI_I;
                    cTRI(Tti)=cTRI_II;

                    indI=[Ttii TRI(j,1) (lUV+i)];
                    indII=[Ttii TRI(j,2) (lUV+i)];
                    
                    TRIcornerI=[TRIcorner(Tti,1),2,TRIcorner(j,2)];
                    TRIcornerII=[TRIcorner(Tti,2),2,TRIcorner(j,3)];
                    
                    TRItricornerI=[TRItricorner(Tti,1),Tti,TRItricorner(j,2)];
                    TRItricornerII=[TRItricorner(Tti,2),j,TRItricorner(j,3)];
                    
                    if TRItricorner(j,3)>0
                        TRItricorner(TRItricorner(j,3),(4-TRIcorner(j,3)))=Tti;
                    end
                    if TRItricorner(Tti,1)>0
                        TRItricorner(TRItricorner(Tti,1),(4-TRIcorner(Tti,1)))=j;
                    end
                    if TRItricorner(j,2)>0
                        TRIcorner(TRItricorner(j,2),(4-TRIcorner(j,2)))=1;
                    end
                    if TRItricorner(Tti,2)>0
                        TRIcorner(TRItricorner(Tti,2),(4-TRIcorner(Tti,2)))=3;
                    end
                    
                    TRI(j,:)=indI;
                    TRI(Tti,:)=indII;
                    
                    TRIcorner(j,:)=TRIcornerI;
                    TRIcorner(Tti,:)=TRIcornerII;
                    
                    TRItricorner(j,:)=TRItricornerI;
                    TRItricorner(Tti,:)=TRItricornerII;                      
                end
            end
        end
    end
end

function comp=compTRI(UV)

AreaTRI=abs(det([(UV(:,2)-UV(:,1)) (UV(:,3)-UV(:,1))]));

les=[norm(UV(:,1)-UV(:,2)),norm(UV(:,1)-UV(:,3)),norm(UV(:,3)-UV(:,2))];

if AreaTRI>1e-12;
    comp=sqrt(max(les)/min(les))*AreaTRI;
else
    comp=0;
end

function uv=triIPoints(UV)

[so,ind]=sort([norm(UV(:,1)-UV(:,2)),norm(UV(:,1)-UV(:,3)),norm(UV(:,3)-UV(:,2))]);

alpha=0.6+0.3*rand;

if ind(1)==1
    uv=alpha*mean(UV(:,[1 2 3 3]),2)+(1-alpha)*mean(UV(:,[1,2]),2);
elseif ind(1)==2
    uv=alpha*mean(UV(:,[1 2 2 3]),2)+(1-alpha)*mean(UV(:,[1,3]),2);
else
    uv=alpha*mean(UV(:,[1 1 2 3]),2)+(1-alpha)*mean(UV(:,[2,3]),2);
end
