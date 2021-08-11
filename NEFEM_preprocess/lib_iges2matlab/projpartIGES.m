function [model,UV,srfind,srfDerivind,srfDer,numpoints]=projpartIGES(ParameterData,R,normal,pdir,sdir,dp,nppos,npneg,nspos,nsneg,plot_flag,fignr)
% PROJPARTIGES returns points of projections on surfaces from an IGES-file.
% 
% Usage:
% 
% [model,UV,srfind,srfDerivind,srfDer,numpoints]=projpartIGES(ParameterData,...
%            R,normal,pdir,sdir,dp,nppos,npneg,nspos,nsneg,plot_flag,fignr)
% 
% Input:
% 
% ParameterData - Parameter data from IGES file. ParameterData
%                 is one of the output from IGES2MATLAB.
% R - The projection origo.
% normal - The projection normal. The direction of normal is toward the surface.
% pdir - The first (primary) direction in which projection points lies.
% sdir - The second (secondary) direction in which projection points lies.
% dp - the distance between the projected points.
% nppos - number of projections in positive primary direction.
% npneg - number of projections in negative primary direction.
% nspos - number of projections in positive secondary direction.
% nsneg - number of projections in negative secondary direction.
% plot_flag - [0/1] =0, do not a plot (default)
%                   =1, do a plot of points of projection and corresponding
%                       surface.
% fignr - Figure number for the plot if plot_flag==1
% 
% Output:
% 
% model - points of projetions.
% UV - the parameter values for corresponding model point from original surface.
% srfind - The index of surface in ParameterData for corresponding model point. srfind(i)==0 means no projection.
% srfDerivind - The index of surface derivatives in srfDer for corresponding model point.
% srfDer - Cell array with surface first and second derivative for all model points.
%               
% m-file can be downloaded for free at
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=
% 13253&objectType=file
% 
% projpartIGES version 1.3
% 
% written by Per Bergström 2007-01-08
%

if nargin<12
    fignr=1;
    if nargin<11
        plot_flag=0;
        if nargin<10
            error('projIGES must have 10 input arguments!');
        end
    end
end

if not(iscell(ParameterData))
    error('ParameterData must be a cell array!');
end

[mR,nR]=size(R);

if mR<nR
    R=R';
    [mR,nR]=size(R);
end

[mnormal,nnormal]=size(normal);

if mnormal<nnormal
    normal=normal';
    [mnormal,nnormal]=size(normal);
end

[mpdir,npdir]=size(pdir);

if mpdir<npdir
    pdir=pdir';
    [mpdir,npdir]=size(pdir);
end

[msdir,nsdir]=size(sdir);

if msdir<nsdir
    sdir=sdir';
    [msdir,nsdir]=size(sdir);
end

if not(and(and(mpdir==3,mnormal==3),and(msdir==3,mR==3)))
    error('Length of R, normal, pdir and sdir must be 3!');
end

if length(dp)~=1
    error('dp must be a scalar!');
elseif dp<eps
    error('dp must be larger than eps!');
end

if or(isempty(nppos),length(nppos)>1)
    nppos=0;
elseif nppos<0
    nppos=0;
else
    nppos=round(nppos);
end

if or(isempty(npneg),length(npneg)>1)
    npneg=0;
elseif npneg<0
    npneg=0;
else
    npneg=round(npneg);
end

if or(isempty(nspos),length(nspos)>1)
    nspos=0;
elseif nspos<0
    nspos=0;
else
    nspos=round(nspos);
end

if or(isempty(nsneg),length(nsneg)>1)
    nsneg=0;
elseif nsneg<0
    nsneg=0;
else
    nsneg=round(nsneg);
end

if isempty(plot_flag)
    plot_flag=0;
elseif not(or(plot_flag==0,plot_flag==1))
    plot_flag=0;
end

if or(isempty(fignr),length(fignr)>1)
    fignr=1;
elseif fignr<1
    fignr=1;
else
    fignr=round(fignr);
end

normal=normal/norm(normal);

pdir=pdir-dot(pdir,normal)*normal;   % primary direction

nopd=norm(pdir);
if nopd<1e-6
    error('pdir can not be parallel to normal');
else
    pdir=pdir/nopd;                  % orthogonal to normal
end

sdirdir=cross(normal,pdir);          
sdirdir=sdirdir/norm(sdirdir);

sdir=dot(sdir,sdirdir)*sdirdir;      % secondary direction

nosd=norm(sdir);
if nosd<1e-6
    error('Illegeal pdir, sdir or normal!');
else
    sdir=sdir/nosd;
end

numpoints=[(nppos+npneg+1),(nspos+nsneg+1)];

pO=[pdir sdir]\(R-pdir*(npneg*dp)-sdir*(nsneg*dp));

[model,UV,srfind,srfDerivind,srfDer,nmodel]=projIGESsub(ParameterData,normal,pdir,sdir,dp,nppos+npneg+1,nspos+nsneg+1,pO);

if plot_flag
    
    figure(fignr),hold on;
    plot3(model(1,srfind>0),model(2,srfind>0),model(3,srfind>0),'b.');
    
    if srfind((nspos+nsneg+1)*npneg+nsneg-1)>0
        plot3(model(1,(nspos+nsneg+1)*npneg+nsneg-1),model(2,(nspos+nsneg+1)*npneg+nsneg-1),model(3,(nspos+nsneg+1)*npneg+nsneg-1),'r.');
    end
    
    clr=[0.53 0.24 0.24];
    lclrf=0.3;
    
    siz=length(ParameterData);
    for i=siz:-1:1
        if ParameterData{i}.type==314
            clr=[ParameterData{i}.cc1 ParameterData{i}.cc2 ParameterData{i}.cc3]/100;
            break
        end
    end
    
    siz2=length(srfDer);
    
    for i=1:siz2
        [P,isSCP,isSup,TRI]=retSrfCrvPnt(1,ParameterData,1,srfDer{i}.supind,100);
        if and(isSCP,not(isSup))
            trisurf(TRI,P(1,:),P(2,:),P(3,:),'FaceColor',clr,'EdgeColor','none','FaceAlpha',0.6);
        end
    end
    
    figure(fignr),hold off;
    
end

function [model,UV,srfind,srfDerivind,srfDer,nmodel]=projIGESsub(ParameterData,normal,pdir,sdir,dp,np1,np2,pO)

nmodel=np1*np2;

model=zeros(3,nmodel);

UV=zeros(2,nmodel);

srfind=zeros(1,nmodel);
srfindsup=zeros(1,nmodel);

normalz=Inf*ones(1,nmodel);

siz=length(ParameterData);

for i=1:siz   % Triangulate each surface and find projection on triangulation

    [PTRI,isSCP,isSup,TRI,UV0,srfind0]=retSrfCrvPnt(1,ParameterData,1,i,500,0);
    
    if and(isSCP,not(isSup))

        mTRI=size(TRI,1);

        for j=1:mTRI

            trip=[pdir sdir]\PTRI(:,TRI(j,:));

            indfloat=(trip-repmat(pO,1,3))./(dp*ones(2,3));

            ind1s=floor(min(indfloat(1,:)))+1;
            ind1e=ceil(max(indfloat(1,:)))+1;
            ind2s=floor(min(indfloat(2,:)))+1;
            ind2e=ceil(max(indfloat(2,:)))+1;
            
            if not(ind1s>np1 | ind2s>np2 | ind1e<1 | ind2e<1)

                if ind1s<1
                    ind1s=1;
                end

                if ind1e>np1
                    ind1e=np1;
                end

                if ind2s<1
                    ind2s=1;
                end

                if ind2e>np2
                    ind2e=np2;
                end

                for ii=ind1s:ind1e
                    for jj=ind2s:ind2e

                        p=pO+[(ii-1);(jj-1)]*dp;

                        rcnd=rcond([(trip(:,1)-trip(:,3)) (trip(:,2)-trip(:,3))]);

                        if rcnd>1e-12
                            ab=[(trip(:,1)-trip(:,3)) (trip(:,2)-trip(:,3))]\(p-trip(:,3));
                            c=1-sum(ab);

                            if min(ab)>=0 & max(ab)<=1 & c>=0

                                Ptemp=ab(1)*PTRI(:,TRI(j,1))+ab(2)*PTRI(:,TRI(j,2))+c*PTRI(:,TRI(j,3));

                                normalztemp=dot(Ptemp,normal);
                                if normalztemp<normalz((ii-1)*np2+jj)
                                    normalz((ii-1)*np2+jj)=normalztemp;
                                    UV(:,(ii-1)*np2+jj)=ab(1)*UV0(:,TRI(j,1))+ab(2)*UV0(:,TRI(j,2))+c*UV0(:,TRI(j,3));
                                    srfind((ii-1)*np2+jj)=srfind0;
                                    srfindsup((ii-1)*np2+jj)=i;
                                    model(:,(ii-1)*np2+jj)=Ptemp;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

[sosrfind,indsoP]=sort(srfind);

sti=1+np1*np2;

%% assume the only surface is type 128

numbder=0;
soi2=0;
for i=1:nmodel
    if and(sosrfind(i)>0,numbder==0)
        sti=i;
        soi=sosrfind(i);
        soi2=soi;
        dnurbs=nrbderivIGES(ParameterData{soi}.nurbs);
        numbder=1;
    elseif sosrfind(i)>soi2
        soi2=sosrfind(i);
        numbder=numbder+1;
    end
end

if numbder>0

    srfDer=cell(1,numbder);
    
    srfDerivind=zeros(1,nmodel);
    
    srfDerind=1;

    if ParameterData{srfind(indsoP(sti))}.type==128
        srfDer{srfDerind}.type=128;
        srfDer{srfDerind}.name='B-NURBS SRF & DERIVATIVE';
        srfDer{srfDerind}.nurbs=ParameterData{soi}.nurbs;
        srfDer{srfDerind}.dnurbs=dnurbs;
        srfDer{srfDerind}.d2nurbs=nrbderiv2IGES(dnurbs);
        srfDer{srfDerind}.supind=srfindsup(indsoP(sti));
    end

end

for i=sti:nmodel % For each model point. Find the projection on the corresponding surface.

    if sosrfind(i)>soi
       
        soi=sosrfind(i);
        
        dnurbs = nrbderivIGES(ParameterData{soi}.nurbs);
        
        srfDerind=srfDerind+1;

        if ParameterData{srfind(indsoP(i))}.type==128
            srfDer{srfDerind}.type=128;
            srfDer{srfDerind}.name='B-NURBS SRF & DERIVATIVE';
            srfDer{srfDerind}.nurbs=ParameterData{soi}.nurbs;
            srfDer{srfDerind}.dnurbs=dnurbs;
            srfDer{srfDerind}.d2nurbs=nrbderiv2IGES(dnurbs);
            srfDer{srfDerind}.supind=srfindsup(indsoP(i));
        end

    end
  
    srfDerivind(indsoP(i))=srfDerind;
    
    if UV(1,indsoP(i))<ParameterData{soi}.u(1)
        UV(1,indsoP(i))=ParameterData{soi}.u(1);
    end
    
    if UV(1,indsoP(i))>ParameterData{soi}.u(2)
        UV(1,indsoP(i))=ParameterData{soi}.u(2);
    end
    
    if UV(2,indsoP(i))<ParameterData{soi}.v(1)
        UV(2,indsoP(i))=ParameterData{soi}.v(1);
    end
    
    if UV(2,indsoP(i))>ParameterData{soi}.v(2)
        UV(2,indsoP(i))=ParameterData{soi}.v(2);
    end     
    
    distmin2p=Inf;
    
    for k=1:50

        [cp,cw] = nrbevalIGES(ParameterData{soi}.nurbs,UV(:,indsoP(i)));

        Psrf = cp./cw;
        
        modmP=model(:,indsoP(i))-Psrf;
        
        distmin2=norm(cross(modmP,normal));
        distmin2ind=1;
        
        if or(distmin2<1e-8,(distmin2p-distmin2)<1e-25)
            break
        end
        
        distmin2p=distmin2;        
 
        [cup,cuw] = nrbevalIGES(dnurbs{1},UV(:,indsoP(i)));
        d1 = (cup-cuw.*Psrf)./cw;
        
        [cvp,cvw] = nrbevalIGES(dnurbs{2},UV(:,indsoP(i)));
        d2 = (cvp-cvw.*Psrf)./cw;      
       
        if rcond([d1 d2 normal])>1e-12

            t=1.2*([d1 d2 normal]\modmP);

            UVt=UV(:,indsoP(i))+[t(1);t(2)];

            if UVt(1)<ParameterData{soi}.u(1)
                UVt=UVt-((UVt(1)-ParameterData{soi}.u(1))/t(1))*[t(1);t(2)];
            end

            if UVt(1)>ParameterData{soi}.u(2)
                UVt=UVt-((UVt(1)-ParameterData{soi}.u(2))/t(1))*[t(1);t(2)];
            end

            if UVt(2)<ParameterData{soi}.v(1)
                UVt=UVt-((UVt(2)-ParameterData{soi}.v(1))/t(2))*[t(1);t(2)];
            end

            if UVt(2)>ParameterData{soi}.v(2)
                UVt=UVt-((UVt(2)-ParameterData{soi}.v(2))/t(2))*[t(1);t(2)];
            end

            numbt=10;
            
            for kk=0:3

                lint=linspace(0,(0.3/(numbt-1))^kk,numbt);
                
                for j=2:numbt

                    UVtemp=(1-lint(j))*UV(:,indsoP(i))+lint(j)*UVt;
                    Ptemp=nrbevalIGES(ParameterData{soi}.nurbs,UVtemp);

                    distmin2temp=norm(cross((model(:,indsoP(i))-Ptemp),normal));

                    if distmin2temp<distmin2
                        distmin2=distmin2temp;
                        distmin2ind=j;
                    end
                end
                
                if distmin2ind>1
                    break
                end

            end

            UV(:,indsoP(i))=(1-lint(distmin2ind))*UV(:,indsoP(i))+lint(distmin2ind)*UVt;
        end
    end

    model(:,indsoP(i))=nrbevalIGES(ParameterData{soi}.nurbs,UV(:,indsoP(i)));
    
end

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

function dnurbs = nrbderivIGES(nurbs)

dnurbs=cell(1,2);

dnurbs{1}.form='B-NURBS';
dnurbs{1}.dim=4;
dnurbs{1}.number=[(nurbs.number(1)-1),nurbs.number(2)];
dnurbs{1}.coefs=(nurbs.order(1)-1)*(nurbs.coefs(:,2:nurbs.number(1),:)-nurbs.coefs(:,1:(nurbs.number(1)-1),:))./repmat((nurbs.knots{1}((1+nurbs.order(1)):(nurbs.number(1)+nurbs.order(1)-1))-nurbs.knots{1}(2:nurbs.number(1))),[4,1,nurbs.number(2)]);
dnurbs{1}.knots={nurbs.knots{1}(2:(end-1)) nurbs.knots{2}};
dnurbs{1}.order=[(nurbs.order(1)-1) nurbs.order(2)];

dnurbs{2}.form='B-NURBS';
dnurbs{2}.dim=4;
dnurbs{2}.number=[nurbs.number(1),(nurbs.number(2)-1)];

dcoefsv=zeros(4,nurbs.number(1),nurbs.number(2)-1);
reshkn2=reshape(nurbs.knots{2}((1+nurbs.order(2)):(nurbs.number(2)+nurbs.order(2)-1))-nurbs.knots{2}(2:nurbs.number(2)),[1,1,(nurbs.number(2)-1)]);
for ii=1:4
    for jj=1:nurbs.number(1)
        dcoefsv(ii,jj,:)=(nurbs.order(2)-1)*(nurbs.coefs(ii,jj,2:nurbs.number(2))-nurbs.coefs(ii,jj,1:(nurbs.number(2)-1)))./reshkn2;
    end
end

dnurbs{2}.coefs=dcoefsv;
dnurbs{2}.knots={nurbs.knots{1} nurbs.knots{2}(2:(end-1))};
dnurbs{2}.order=[nurbs.order(1) (nurbs.order(2)-1)];

function d2nurbs = nrbderiv2IGES(dnurbs)

warning off;   % srf perhaps not diffrentiable 2 times

d2nurbs11_12 = nrbderivIGES(dnurbs{1});
d2nurbs21_22 = nrbderivIGES(dnurbs{2});

warning on;

d2nurbs=cell(1,3);

d2nurbs{1}=d2nurbs11_12{1};
d2nurbs{2}=d2nurbs11_12{2};
d2nurbs{3}=d2nurbs21_22{2};

%reduce roundoff error

d2nurbs{2}.coefs=0.5*(d2nurbs11_12{2}.coefs+d2nurbs21_22{1}.coefs);