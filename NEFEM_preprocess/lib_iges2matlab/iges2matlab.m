function [ParameterData,EntityType,numEntityType,unknownEntityType,numunknownEntityType]=iges2matlab(igsfile)
% IGES2MATLAB converts the parameter data in the IGES-file igsfile
% to MATLAB format.
% 
% Usage:
% 
% [ParameterData,EntityType,numEntityType,unknownEntityType,numunknownEntityType]=iges2matlab(igsfile)
% 
% Input:
% 
% igsfile - IGES file
% 
% Output:
% 
% ParameterData - cell array with Parameter Data from igsfile
% EntityType - vector with entities in igsfile converted to matlab
% numEntityType - vector with number of entities belonging to EntityType
% unknownEntityType - vector with unknown entities for iges2matlab
% numunknownEntityType - vector with number of unknown entities
%                        belonging to unknownEntityType
% 
% For entity type 126 and 128 ParameterData also contains a nurbs
% representation same as in NURBS toolbox. NURBS toolbox can be downloaded
% for free at
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=
% 312&objectType=file
% 
% ParameterData also contains other useful information for usage in other
% functions. For curves the length is given as a parameter in ParameterData.
% superior is another parameter for curves and surfaces. For curves superior=1
% means that they are defined in the parameter space for a surface. superior=0
% means that they are defined in the 3D-space. For surfaces superior=1 means
% that their domain is limited by one or more closed curves. superior=0 means
% that their domain is the domain given in the ParameterData.
% 
% For other parameters see the IGES specificaton version5x.pdf found at
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=
% 12087&objectType=file
% 
% All pointers in ParameterData points to the index 
% in cell array ParameterData.
% 
% This version can not handle all possible IGES entities.
% 
% Example:
% 
% [ParameterData,EntityType,numEntityType,unknownEntityType]=iges2matlab('tool.igs');
% 
% will convert the parameter data in tool.igs to MATLAB format.
% 
% m-file can be downloaded for free at
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=
% 13253&objectType=file
% 
% iges2matlab version 1.3
% 
% written by Per Bergstr�m 2007-01-08
%

[fid,msg]=fopen(igsfile);

if fid==-1
    error(msg);
end

c = fread(fid,'uint8=>uint8')';

fclose(fid);

nwro=sum((c((81:82))==10))+sum((c((81:82))==13));

edfi=nwro-sum(c(((end-1):end))==10)-sum(c(((end-1):end))==13);

siz=length(c);

ro=round((siz+edfi)/(80+nwro));

if rem((siz+edfi),(80+nwro))~=0
    error('Input file must be an IGES-file!');
end

roind=1:ro;

SGDPT=c(roind*(80+nwro)-7-nwro);

Sfind=SGDPT==83;
Gfind=SGDPT==71;
Dfind=SGDPT==68;
Pfind=SGDPT==80;
Tfind=SGDPT==84;

sumSfind=sum(Sfind);
sumGfind=sum(Gfind);
sumDfind=sum(Dfind);
sumPfind=sum(Pfind);
sumTfind=sum(Tfind);

%------S Line information (The initial line to get things started)---------
for i=roind(Sfind)
    disp(char(c(((i-1)*(80+nwro)+1):(i*(80+nwro)-8-nwro))));
end

%---------------G Line information  (Header infomation)--------------------

G=cell(1,25);
Gstr=zeros(1,72*sumGfind);
j=1;
for i=roind(Gfind)
    Gstr(((j-1)*72+1):(j*72))=c(((i-1)*(80+nwro)+1):(i*(80+nwro)-8-nwro));
    j=j+1;
end

if and(Gstr(1)==49,Gstr(2)==72)
    G{1}=Gstr(3);
    st=4;
else
    G{1}=44;
    st=1;
end

if and(Gstr(st+1)==49,Gstr(st+2)==72)
    G{2}=Gstr(st+3);
    st=st+4;
else
    G{2}=59;
    st=st+1;
end

le=length(Gstr);

for i=3:25
    for j=(st+1):le
        if or(Gstr(j)==G{1},Gstr(j)==G{2})
            break
        end
    end
    G{i}=Gstr((st+1):(j-1));
    st=j;
end

for i=[3 4 5 6 12 15 18 21 22 25]   %string
    stind=1;
    for j=1:length(G{i})
        if G{i}(j)~=32
            stind=j;
            break
        end 
    end
    
    for j=stind:length(G{i})
        if G{i}(j)==72
            stind=j+1;
            break
        end 
    end        
    
    endind=length(G{i});
    for j=length(G{i}):-1:1
        if G{i}(j)~=32
            endind=j;
            break
        end 
    end
    G{i}=G{i}(stind:endind);    
end

for i=[7 8 9 10 11 13 14 16 17 19 20 23 24]   %num
    G{i}=str2num(char(G{i}));
end

%--D Line information (Data information) & P Line information (All data)---

noent=round(sumDfind/2);

ParameterData=cell(1,noent);

roP=sumSfind+sumGfind+sumDfind;

entty=zeros(1,520);
entunk=zeros(1,520);

entiall=0;

for i=(sumSfind+sumGfind+1):2:(sumSfind+sumGfind+sumDfind-1)
    
    entiall=entiall+1;
    
    Dstr=c(((i-1)*(80+nwro)+1):(i*(80+nwro)-8-nwro));
    
    type=str2num(char(Dstr(1:8)));
    
    ParameterData{entiall}.type=type;
    
    Pstart=str2num(char(Dstr(9:16)))+roP;
    
    if i==roP-1
        Pend=ro-sumTfind;
    else
        Pend=str2num(char(c(((i+1)*(80+nwro)+9):((i+1)*(80+nwro)+16))))+roP-1;
    end
    
    Pstr=zeros(1,64*(Pend-Pstart+1));
    j=1;
    for k=Pstart:Pend
        Pstr(((j-1)*64+1):(j*64))=c(((k-1)*(80+nwro)+1):(k*(80+nwro)-16-nwro));
        j=j+1;
    end
    
    Pstr(Pstr==G{1})=44;
    Pstr(Pstr==G{2})=59;
    
    Pvec=str2num(char(Pstr));
    
    % SURFACES
    
    if type==128
        
        entty(type)=entty(type)+1;
        
        A=1+Pvec(2)+Pvec(4);
        B=1+Pvec(3)+Pvec(5);
        C=(Pvec(2)+1)*(Pvec(3)+1);
        
        ParameterData{entiall}.name='B-NURBS SRF';
        
        ParameterData{entiall}.superior=0;
        
        ParameterData{entiall}.k1=Pvec(2);
        ParameterData{entiall}.k2=Pvec(3);
        ParameterData{entiall}.m1=Pvec(4);
        ParameterData{entiall}.m2=Pvec(5);
        
        ParameterData{entiall}.prop1=Pvec(6);
        ParameterData{entiall}.prop2=Pvec(7);
        ParameterData{entiall}.prop3=Pvec(8);
        ParameterData{entiall}.prop4=Pvec(10);
        ParameterData{entiall}.prop5=Pvec(11);
        
        ParameterData{entiall}.s=Pvec(11:(11+A));
        ParameterData{entiall}.t=Pvec((12+A):(12+A+B));
        
        ParameterData{entiall}.w=reshape(Pvec((13+A+B):(12+A+B+C)),Pvec(2)+1,Pvec(3)+1);
        
        ParameterData{entiall}.p=zeros(3,Pvec(2)+1,Pvec(3)+1);
        ParameterData{entiall}.p=reshape(Pvec((13+A+B+C):(12+A+B+4*C)),3,Pvec(2)+1,Pvec(3)+1);
        
        ParameterData{entiall}.u=zeros(1,2);
        ParameterData{entiall}.u(1)=Pvec(13+A+B+4*C);
        ParameterData{entiall}.u(2)=Pvec(14+A+B+4*C);
        
        ParameterData{entiall}.v=zeros(1,2);
        ParameterData{entiall}.v(1)=Pvec(15+A+B+4*C);
        ParameterData{entiall}.v(2)=Pvec(16+A+B+4*C);
        
        % NURBS surface
        
        ParameterData{entiall}.nurbs.form='B-NURBS';
        
        ParameterData{entiall}.nurbs.dim=4;
        
        ParameterData{entiall}.nurbs.number=zeros(1,2);
        ParameterData{entiall}.nurbs.number(1)=Pvec(2)+1;
        ParameterData{entiall}.nurbs.number(2)=Pvec(3)+1;
        
        ParameterData{entiall}.nurbs.coefs=zeros(4,Pvec(2)+1,Pvec(3)+1);
        ParameterData{entiall}.nurbs.coefs(4,:,:)=reshape(Pvec((13+A+B):(12+A+B+C)),Pvec(2)+1,Pvec(3)+1);
        ParameterData{entiall}.nurbs.coefs(1:3,:,:)=reshape(Pvec((13+A+B+C):(12+A+B+4*C)),3,Pvec(2)+1,Pvec(3)+1);
        
        ParameterData{entiall}.nurbs.coefs(1:3,:,:)=repmat(ParameterData{entiall}.nurbs.coefs(4,:,:),3,1).*ParameterData{entiall}.nurbs.coefs(1:3,:,:);
        
        ParameterData{entiall}.nurbs.knots=cell(1,2);
        ParameterData{entiall}.nurbs.knots{1}=Pvec(11:(11+A));
        ParameterData{entiall}.nurbs.knots{2}=Pvec((12+A):(12+A+B));
        
        ParameterData{entiall}.nurbs.order=zeros(1,2);
        ParameterData{entiall}.nurbs.order(1)=Pvec(4)+1;
        ParameterData{entiall}.nurbs.order(2)=Pvec(5)+1;
        
    elseif type==144
        
        entty(type)=entty(type)+1;
        
        ParameterData{entiall}.name='TRIMMED SURFACE';
        
        ParameterData{entiall}.pts=round((Pvec(2)+1)/2);
        
        ParameterData{entiall}.n1=Pvec(3);
        ParameterData{entiall}.n2=Pvec(4);
        
        if Pvec(5)~=0
            ParameterData{entiall}.pto=round((Pvec(5)+1)/2);
        else
            ParameterData{entiall}.pto=0;
        end
        
        ParameterData{entiall}.pti=round((Pvec(6:(5+Pvec(4)))+1)/2);
        
        % CURVES
        
    elseif type==126
        
        entty(type)=entty(type)+1;
        
        N=1+Pvec(2)-Pvec(3);
        A=1+Pvec(2)+Pvec(3);
        
        ParameterData{entiall}.name='B-NURBS CRV';
        
        ParameterData{entiall}.superior=0;
        
        ParameterData{entiall}.k=Pvec(2);
        ParameterData{entiall}.m=Pvec(3);
        
        ParameterData{entiall}.prop1=Pvec(4);
        ParameterData{entiall}.prop2=Pvec(5);
        ParameterData{entiall}.prop3=Pvec(6);
        ParameterData{entiall}.prop4=Pvec(7);
        
        ParameterData{entiall}.t=Pvec(8:(8+A));
        
        ParameterData{entiall}.w=Pvec((9+A):(9+A+Pvec(2)));
        
        ParameterData{entiall}.p=reshape(Pvec((10+A+Pvec(2)):(12+A+4*Pvec(2))),3,Pvec(2)+1);
        
        
        ParameterData{entiall}.v=zeros(1,2);
        ParameterData{entiall}.v(1)=Pvec(13+A+4*Pvec(2));
        ParameterData{entiall}.v(2)=Pvec(14+A+4*Pvec(2));
        
        if Pvec(4)
            ParameterData{entiall}.xnorm=Pvec(15+A+4*Pvec(2));
            ParameterData{entiall}.ynorm=Pvec(16+A+4*Pvec(2));
            ParameterData{entiall}.znorm=Pvec(17+A+4*Pvec(2));
        else
            ParameterData{entiall}.xnorm=0;
            ParameterData{entiall}.ynorm=0;
            ParameterData{entiall}.znorm=0;
        end
        
        % NURBS curve
        
        ParameterData{entiall}.nurbs.form='B-NURBS';
        
        ParameterData{entiall}.nurbs.dim=4;
        
        ParameterData{entiall}.nurbs.number=Pvec(2)+1;
        
        ParameterData{entiall}.nurbs.coefs=zeros(4,Pvec(2)+1);
        ParameterData{entiall}.nurbs.coefs(4,:)=Pvec((9+A):(9+A+Pvec(2)));
        
        ParameterData{entiall}.nurbs.coefs(1:3,:)=reshape(Pvec((10+A+Pvec(2)):(12+A+4*Pvec(2))),3,Pvec(2)+1);
        
        ParameterData{entiall}.nurbs.coefs(1:3,:)=repmat(ParameterData{entiall}.nurbs.coefs(4,:),3,1).*ParameterData{entiall}.nurbs.coefs(1:3,:);
        
        ParameterData{entiall}.nurbs.order=Pvec(3)+1;
        
        ParameterData{entiall}.nurbs.knots=Pvec(8:(8+A));
        
        nup=500;
        p = nrbevalIGES(ParameterData{entiall}.nurbs,linspace(ParameterData{entiall}.v(1),ParameterData{entiall}.v(2),nup));
        ParameterData{entiall}.length=sum(sqrt(sum((p(:,1:(nup-1))-p(:,2:nup)).^2,1)));
        
        
    elseif type==110
        
        entty(type)=entty(type)+1;
        
        ParameterData{entiall}.name='LINE';
        
        ParameterData{entiall}.superior=0;
        
        ParameterData{entiall}.form=str2num(char(c((i*(80+nwro)+65):(i*(80+nwro)+72))));
        
        ParameterData{entiall}.p1=Pvec(2:4)';
        ParameterData{entiall}.x1=Pvec(2);
        ParameterData{entiall}.y1=Pvec(3);
        ParameterData{entiall}.z1=Pvec(4);
        
        ParameterData{entiall}.p2=Pvec(5:7)';
        ParameterData{entiall}.x2=Pvec(5);
        ParameterData{entiall}.y2=Pvec(6);
        ParameterData{entiall}.z2=Pvec(7);
        
        ParameterData{entiall}.length=norm(Pvec(2:4)-Pvec(5:7));
        
    elseif type==102
        
        entty(type)=entty(type)+1;
        
        ParameterData{entiall}.name='COMPOSITE';
        
        ParameterData{entiall}.n=Pvec(2);
        ParameterData{entiall}.de=round((Pvec(3:(2+Pvec(2)))+1)/2);
        
        ParameterData{entiall}.lengthcnt=zeros(1,Pvec(2));
        ParameterData{entiall}.length=0;
        
    elseif type==142
        
        entty(type)=entty(type)+1;
        
        ParameterData{entiall}.name='CRV ON A PARAMETRIC SURFACE';
        
        ParameterData{entiall}.crtn=Pvec(2);
        ParameterData{entiall}.sptr=round((Pvec(3)+1)/2);
        ParameterData{entiall}.bptr=round((Pvec(4)+1)/2);
        ParameterData{entiall}.cptr=round((Pvec(5)+1)/2);
        ParameterData{entiall}.pref=Pvec(6);
        
    elseif type==116
        
        entty(type)=entty(type)+1;
        
        ParameterData{entiall}.name='POINT';
        
        ParameterData{entiall}.p=Pvec(2:4)';
        
        ParameterData{entiall}.x=Pvec(2);
        ParameterData{entiall}.y=Pvec(3);
        ParameterData{entiall}.z=Pvec(4);
        
        %%%%%%%%%%%%%%%%%ParameterData{entiall}.ptr=round((Pvec(5)+1)/2);
        
        % OTHER
        
    elseif type==314
        
        entty(type)=entty(type)+1;
        
        ParameterData{entiall}.name='COLOR';
        
        inn=find(or(Pstr==44,Pstr==59));
        
        ParameterData{entiall}.cc1=str2num(char(Pstr((inn(1)+1):(inn(2)-1))));
        ParameterData{entiall}.cc2=str2num(char(Pstr((inn(2)+1):(inn(3)-1))));
        ParameterData{entiall}.cc3=str2num(char(Pstr((inn(3)+1):(inn(4)-1))));
        
        if length(inn)>4
            inn2=find(Pstr(1:(inn(5)-1))==72);
            if isempty(inn2)
                ParameterData{entiall}.cname='';
            else
                ParameterData{entiall}.cname=char(Pstr((inn2(1)+1):(inn(5)-1)));
            end
        else
            ParameterData{entiall}.cname='';
        end
        
    else
        
        ParameterData{entiall}.name='Unknown type!';
        entunk(type)=entunk(type)+1;
        
    end
end

ent_ind=1:520;

EntityType=ent_ind(entty>0);

numEntityType=entty(entty>0);

unknownEntityType=ent_ind(entunk>0);

numunknownEntityType=entunk(entunk>0);


% Info to plotIGES

for i=1:noent
    if ParameterData{i}.type==144
        ParameterData{ParameterData{i}.pts}.superior=1;
        
        if ParameterData{i}.n1
            ParameterData=dosuperior(ParameterData,ParameterData{i}.pto);
        end
        
        for j=1:ParameterData{i}.n2
            ParameterData=dosuperior(ParameterData,ParameterData{i}.pti(j));
        end
    elseif ParameterData{i}.type==102
        
        for j=1:ParameterData{i}.n
            ParameterData{i}.lengthcnt(j)=ParameterData{ParameterData{i}.de(j)}.length;
        end
        
        ParameterData{i}.length=sum(ParameterData{i}.lengthcnt);
    end
end

% Recursive define function

function ParameterData=dosuperior(ParameterData,ii)

ty=ParameterData{ii}.type;

if ty==126
    ParameterData{ii}.superior=1;
elseif ty==110
    ParameterData{ii}.superior=1;
elseif ty==102
    for k=1:ParameterData{ii}.n
        ParameterData=dosuperior(ParameterData,ParameterData{ii}.de(k));
    end
elseif ty==142
    % only bptr, not cptr
    ParameterData=dosuperior(ParameterData,ParameterData{ii}.bptr);
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






        