function data = run_simulation_DG(data,handle)

U = data.U;
mesh = data.mesh;
T = mesh.T;
BC = data.BC;
referenceElement = mesh.referenceElement;
elementalMatricesInfo = data.elementalMatricesInfo;
t = data.t;
infoFaces = data.infoFaces;
incidentWave = data.ip;
optStore= data.sp.optStore;
tolSteadyState= data.sp.sstol;
nOfCyclesMax = data.sp.nOfCyclesMax;
dt = data.sp.dt_stab/data.sp.stabParam;
if any(strcmp(data.PML(5,:),'on'))
    PML = data.PMLabsorptionValue;
else
    PML = [];
end
c = data.bottom.c;
cg = data.bottom.cg;
ccg = c.*cg;
% Information for the time marching process--------------------------------
tEnd = incidentWave.period*nOfCyclesMax;
% nStep = round(tEnd/dt);

% Control the time-harmonic steady state-----------------------------------
boolSteadyState = 0;
errSteadyState = zeros(1,nOfCyclesMax);
% reference solution initialization----------------------------------------
referenceSolution = zeros(2); %needed to run Is... .cxx
% Cycle counter -----------------------------------------------------------
nStepByCycle = round(tEnd/(dt*nOfCyclesMax));
nStep = nStepByCycle*nOfCyclesMax;
iCycle = int32(0);

sstBoundary = data.mesh.boundaryNames{data.sp.sstBoundary};
% figure(2)
% postprocessSurfaceFEM(mesh, referenceElement, U, 4);

dt = tEnd/nStep; % TIME STEP!
nOfPlots = max(c)*dt*nStep/data.sp.spaceXtimeFrame;
nDib = round(nStep/nOfPlots);
if nDib == 0
    nDib = 1;
end
setOutput({'Starting time computation'},handle)
cpuTimeStart = cputime;
%% Loop for time steps
cpuT = cputime;
for iStep = 1:1
    % RK4------------------------------------------------------------------
    F1 = hyperbolicF_FEM(U,mesh,referenceElement,elementalMatricesInfo,...
        t,infoFaces,incidentWave,BC,c,ccg,PML);
    Utmp = U - (dt/2)*F1;
    t = t + dt/2;
    F2 = hyperbolicF_FEM(Utmp,mesh,referenceElement,elementalMatricesInfo,...
        t,infoFaces,incidentWave,BC,c,ccg,PML);
    Utmp = U - (dt/2)*F2;
    F3 = hyperbolicF_FEM(Utmp,mesh,referenceElement,elementalMatricesInfo,...
        t,infoFaces,incidentWave,BC,c,ccg,PML);
    Utmp = U - dt*F3;
    t = iStep*dt;
    F4 = hyperbolicF_FEM(Utmp,mesh,referenceElement,elementalMatricesInfo,...
        t,infoFaces,incidentWave,BC,c,ccg,PML);
    %    U = U - (dt/6)*(F1+2*F2+2*F3+F4);
    for iElem = 1:size(U,3)
        Te = T(iElem,:);
        c_el = c(Te);
        cg_el = cg(Te);
        U(:,1,iElem) = U(:,1,iElem) - (dt/6)*c_el./cg_el.*(F1(:,1,iElem)+2*F2(:,1,iElem)+2*F3(:,1,iElem)+F4(:,1,iElem));
    end
    U(:,2:end,:) = U(:,2:end,:) - (dt/6)*(F1(:,2:end,:)+2*F2(:,2:end,:)+2*F3(:,2:end,:)+F4(:,2:end,:));

    if(mod(iStep,nDib)==0 || iStep==nStep)
        disp(sprintf('t=%1.3f \t CPU Time: %3.3f', t,cputime-cpuT));
        cpuT = cputime;
        if  ~isempty(find(abs(U)>1000,1))
            disp('Divergence detected');
            break;
        end
        storeSolution(U,incidentWave,t,mesh,referenceElement,optStore);
    end
    % Control the time-harmonic steady state-------------------------------
    if  mod(iStep,nStepByCycle)==0
        [iCycle,errSteadyState,boolSteadyState,referenceSolution] = IsSteadyStateAchieved...
            (iCycle,errSteadyState, referenceElement,...
            nOfCyclesMax,U,referenceSolution,mesh,tolSteadyState,infoFaces,sstBoundary);
    end
    if boolSteadyState
%         figure(2)
%         clf
%         postprocessSurfaceFEM(mesh, referenceElement, U, 4);
%         caxis([-1 1])
%         pause(0.1)
        break;
    end
end
setOutput({'Done'},handle)
data.cputime.run_time = cputime - cpuTimeStart;
errSteadyState = errSteadyState(1:iCycle);
data.solution = reshape(U(:,4,:),numel(T),1);
data.solutionData.errorSteadyState = errSteadyState;
data.solutionData.finalTime = t;
data.solutionData.cicles2converge = iCycle;
data.solutionData.dt = dt;


%% FEM
function F = hyperbolicF_FEM(U,mesh,referenceElement,elementalMatricesInfo,...
    t,infoFaces,incidentWave,BC,c,ccg,PML)

F = berkoffFacesTerm(U,mesh,referenceElement,infoFaces,...
    elementalMatricesInfo,t,incidentWave,BC,c,ccg);

F = berkoffDivergenceTerm(F,U,mesh,...
    elementalMatricesInfo,ccg);

F = berkoffForceTerm(F,U,mesh,...
    elementalMatricesInfo,t,incidentWave,PML);

%%
function storeSolution(U,incidentWave,t,mesh,referenceElement,optStore)

% 0=no plots (Clonetroop)
% 1=plot reflected potential
% 2=plot reflected + incident potential
% 3=store reflected potential

if optStore==1
    figure(2)
    clf
    postprocessSurfaceFEM(mesh, referenceElement, U, 4);
    caxis([-1 1])
    pause(0.1)
end

if optStore==2
    figure(2)
    clf
    phi0 = incidentPotential(mesh,incidentWave,t) ;
    postprocessSurfaceFEM(mesh, referenceElement, phi0 + U(:,4,:));
    caxis([-1 1])
    pause(0.1)

end


