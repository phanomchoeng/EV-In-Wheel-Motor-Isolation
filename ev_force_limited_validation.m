clear; close all; clc;

% REVIEWER1_FORCE_LIMITED_VALIDATION
% Standalone reproduction of the additional analyses requested by Reviewer 1:
% 1) dynamic load coefficient (DLC),
% 2) actuator saturation occurrence and mechanical power, and
% 3) force-limited speed sweep for the revised manuscript Figure 10.
%
% The script does not require MAT files produced by other scripts. It creates
% reviewer1_data and reviewer1_figures beside this file. MATLAB R2025b was
% used for the reported results. Control System Toolbox is required.

rootDir = fileparts(mfilename('fullpath'));
dataDir = fullfile(rootDir,'reviewer1_data');
figDir = fullfile(rootDir,'reviewer1_figures');
if ~exist(dataDir,'dir'), mkdir(dataDir); end
if ~exist(figDir,'dir'), mkdir(figDir); end

% Mass-matched configurations and selected passive design reported in the
% manuscript. The total vehicle corner mass is 400 kg in every case.
pCentral = struct('ms',370,'mu',30,'kt',192000,'bs',1400,'keq',20400);
pRigid = struct('ms',330,'mu',70,'kt',192000,'bs',1400,'keq',20400);
pIso = struct('ms',330,'mi',40,'mu',30,'kt',192000,'bs',1400, ...
    'ku',38118.90,'kl',43886.78,'ci',1774.35);

% Full-state LQR gains reported in Section 3.2. State order:
% x = [z_u, dz_u, z_i, dz_i, z_s, dz_s]^T.
Kcomfort = [-337.62,-1344.14,-36638.44,-1740.99,33580.61,1591.96];
Kbalanced = [-1374.98,-1156.20,-31238.44,-1580.14,28206.58,1092.89];
controllerNames = ["Comfort LQR";"Balanced LQR"];
gains = {Kcomfort,Kbalanced};
selectedLimits = [1500;1000];

set(groot,'defaultAxesFontName','Times New Roman', ...
    'defaultAxesFontSize',10,'defaultAxesLineWidth',0.8, ...
    'defaultLineLineWidth',1.35,'defaultFigureColor','w');

Mcorner = pIso.ms + pIso.mi + pIso.mu;
g = 9.81;
Fz0 = Mcorner*g;

[Ac,Brc] = conventional_model(pCentral);
[Ar,Brr] = conventional_model(pRigid);
[Ai,Bri,BFi] = intermediate_model(pIso);
sysCentral = passive_output_system(Ac,Brc,pCentral,false);
sysRigid = passive_output_system(Ar,Brr,pRigid,false);
sysIso = passive_output_system(Ai,Bri,pIso,true);

%% Thirty-road ensemble and dynamic load coefficient (Table 5)
dtR = 0.002; tR = (0:dtR:60).'; roadEval = tR>=5; nRoad = 30;
ens = nan(nRoad,5,3); % RMS acceleration, travel (mm), tire force
for s=1:nRoad
    zr = iso8608_profile(tR,20,20261000+s);
    yc = lsim(sysCentral,zr,tR);
    yr = lsim(sysRigid,zr,tR);
    yi = lsim(sysIso,zr,tR);
    yComfort = simulate_saturated_diagnostics( ...
        Ai,Bri,BFi,Kcomfort,pIso,zr,tR,selectedLimits(1));
    yBalanced = simulate_saturated_diagnostics( ...
        Ai,Bri,BFi,Kbalanced,pIso,zr,tR,selectedLimits(2));
    YY = {yc,yr,yi,yComfort,yBalanced};
    for j=1:5
        m = response_metrics(YY{j}(roadEval,:));
        ens(s,j,:) = [m.rmsAcc,1e3*m.rmsTravel,m.rmsTire];
    end
end
ensDLC = ens(:,:,3)/Fz0;
configNames = ["Centralized drive";"Rigid wheel-motor"; ...
    "Passive isolated motor";"Comfort LQR, +/-1500 N"; ...
    "Balanced LQR, +/-1000 N"];
dlcStats = nan(5,4);
for j=1:5
    v = ensDLC(:,j);
    dlcStats(j,:) = [mean(v),std(v),prctile(v,5),prctile(v,95)];
end
TDLC = table(configNames,dlcStats(:,1),dlcStats(:,2),dlcStats(:,3),dlcStats(:,4), ...
    'VariableNames',{'Configuration','Mean_DLC','Std_DLC','P05_DLC','P95_DLC'});
writetable(TDLC,fullfile(dataDir,'Reviewer1_Table5_DLC.csv'));

%% Selected force limits: bump saturation and power
dtB = 0.001; tB = (0:dtB:10).'; t0 = 1; hb = 0.05; Tb = 0.20;
zrB = bump_input(tB,t0,hb,Tb);
bumpRows = nan(2,8);
for j=1:2
    [y,aux] = simulate_saturated_diagnostics(Ai,Bri,BFi,gains{j},pIso,zrB,tB,selectedLimits(j));
    id = tB>=t0;
    m = response_metrics(y(id,:));
    d = actuator_diagnostics(aux,id,dtB);
    bumpRows(j,:) = [m.rmsAcc,m.peakAcc,m.maxTravel,m.peakForce, ...
        d.saturationPercent,d.saturationDuration,d.rmsPower,d.peakAbsPower];
end
TBump = table(controllerNames,selectedLimits,bumpRows(:,1),bumpRows(:,2), ...
    1e3*bumpRows(:,3),bumpRows(:,4),bumpRows(:,5),bumpRows(:,6), ...
    bumpRows(:,7),bumpRows(:,8), ...
    'VariableNames',{'Controller','Force_limit_N','RMS_acc_mps2', ...
    'Peak_acc_mps2','Max_travel_mm','Peak_force_N', ...
    'Saturation_percent','Saturation_duration_s','RMS_mech_power_W', ...
    'Peak_abs_mech_power_W'});
writetable(TBump,fullfile(dataDir,'Reviewer1_Selected_Limits_Bump.csv'));

%% Selected force limits: thirty-road saturation and power ensemble
roadRaw = nan(nRoad,2,7);
% Metrics: saturation %, saturation duration, RMS power, peak |power|,
% RMS force, RMS acceleration, DLC.
for s=1:nRoad
    zr = iso8608_profile(tR,20,20261000+s);
    for j=1:2
        [y,aux] = simulate_saturated_diagnostics(Ai,Bri,BFi,gains{j},pIso,zr,tR,selectedLimits(j));
        m = response_metrics(y(roadEval,:));
        d = actuator_diagnostics(aux,roadEval,dtR);
        roadRaw(s,j,:) = [d.saturationPercent,d.saturationDuration, ...
            d.rmsPower,d.peakAbsPower,m.rmsForce,m.rmsAcc,m.rmsTire/Fz0];
    end
    fprintf('Reviewer 1 road diagnostics %d/%d\n',s,nRoad);
end

metricNames = ["Saturation_percent";"Saturation_duration_s"; ...
    "RMS_mech_power_W";"Peak_abs_mech_power_W";"RMS_force_N"; ...
    "RMS_acc_mps2";"DLC"];
summaryRows = nan(2*numel(metricNames),4);
summaryController = strings(2*numel(metricNames),1);
summaryMetric = strings(2*numel(metricNames),1);
rr = 0;
for j=1:2
    for q=1:numel(metricNames)
        rr=rr+1;
        v=squeeze(roadRaw(:,j,q));
        summaryController(rr)=controllerNames(j);
        summaryMetric(rr)=metricNames(q);
        summaryRows(rr,:)=[mean(v),std(v),prctile(v,5),prctile(v,95)];
    end
end
TRoad = table(summaryController,summaryMetric,summaryRows(:,1),summaryRows(:,2), ...
    summaryRows(:,3),summaryRows(:,4), ...
    'VariableNames',{'Controller','Metric','Mean','Standard_deviation','P05','P95'});
writetable(TRoad,fullfile(dataDir,'Reviewer1_Selected_Limits_Road_Ensemble.csv'));

%% Revised manuscript Figure 10: force-limited speed sweep
speedPlot = (10:1:30).';
speedMetric = nan(numel(speedPlot),5,4); % acceleration, travel, tire force, DLC
speedNames = {'Centralized drive','Rigid wheel-motor','Passive isolated motor', ...
    'Comfort LQR, +/-1500 N','Balanced LQR, +/-1000 N'};

for iv=1:numel(speedPlot)
    v=speedPlot(iv);
    % Resetting the same seed preserves the spatial phase realization; only
    % the time mapping changes with vehicle speed.
    zr=iso8608_profile(tR,v,20260925);
    yc=lsim(sysCentral,zr,tR);
    yr=lsim(sysRigid,zr,tR);
    yi=lsim(sysIso,zr,tR);
    yComfort=simulate_saturated_diagnostics(Ai,Bri,BFi,Kcomfort,pIso,zr,tR,1500);
    yBalanced=simulate_saturated_diagnostics(Ai,Bri,BFi,Kbalanced,pIso,zr,tR,1000);
    YY={yc,yr,yi,yComfort,yBalanced};
    for j=1:5
        m=response_metrics(YY{j}(roadEval,:));
        speedMetric(iv,j,:)=[m.rmsAcc,m.rmsTravel,m.rmsTire,m.rmsTire/Fz0];
    end
    fprintf('Reviewer 1 force-limited speed sweep %d/%d (%g m/s)\n',iv,numel(speedPlot),v);
end

varNames={'Vehicle_speed_mps'};
for j=1:5
    stem=matlab.lang.makeValidName(speedNames{j});
    varNames=[varNames,{[stem '_RMS_acc_mps2'],[stem '_RMS_travel_m'], ...
        [stem '_RMS_tire_N'],[stem '_DLC']}]; %#ok<AGROW>
end
speedMatrix=speedPlot;
for j=1:5
    speedMatrix=[speedMatrix,squeeze(speedMetric(:,j,:))]; %#ok<AGROW>
end
TSpeed=array2table(speedMatrix,'VariableNames',varNames);
writetable(TSpeed,fullfile(dataDir,'Reviewer1_Figure10_Force_Limited_Speed.csv'));

colors={[0 0 0],[0.45 0.45 0.45],[0 0.4470 0.7410], ...
    [0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560]};
styles={'--',':','-','--','-.'};
f=figure('Position',[80 70 1280 650]);
tl=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
titles={'RMS acceleration (m s^{-2})','RMS suspension travel (m)','RMS tire force (N)'};
for q=1:3
    a=nexttile(tl); hold(a,'on'); grid(a,'on');
    for j=1:5
        plot(a,speedPlot,speedMetric(:,j,q),'Color',colors{j}, ...
            'LineStyle',styles{j},'DisplayName',speedNames{j});
    end
    xlabel(a,'Vehicle speed (m s^{-1})'); ylabel(a,titles{q}); title(a,titles{q});
    if q==1, legend(a,'Location','best','FontSize',8); end
end
title(tl,'Class C Random-Road Performance with Selected Actuator Force Limits', ...
    'FontWeight','bold');
exportgraphics(tl,fullfile(figDir,'Figure10_Force_Limited_Speed.png'),'Resolution',300);
exportgraphics(tl,fullfile(figDir,'Figure10_Force_Limited_Speed.pdf'),'ContentType','vector');
savefig(f,fullfile(figDir,'Figure10_Force_Limited_Speed.fig'));
close(f);

%% Human-readable summary and MAT archive
save(fullfile(dataDir,'reviewer1_results.mat'),'TDLC','TBump','TRoad', ...
    'TSpeed','roadRaw','speedMetric','Fz0','selectedLimits');
fid=fopen(fullfile(dataDir,'reviewer1_summary.txt'),'w');
fprintf(fid,'Reviewer 1 additional numerical results\n');
fprintf(fid,'Static corner load Fz0 = %.3f N\n\n',Fz0);
fprintf(fid,'Selected-limit bump results\n');
for j=1:2
    fprintf(fid,'%s, limit +/-%.0f N: saturation %.4f%% (%.4f s), RMS power %.3f W, peak |power| %.3f W, peak force %.3f N\n', ...
        controllerNames(j),selectedLimits(j),TBump.Saturation_percent(j), ...
        TBump.Saturation_duration_s(j),TBump.RMS_mech_power_W(j), ...
        TBump.Peak_abs_mech_power_W(j),TBump.Peak_force_N(j));
end
fprintf(fid,'\nThirty-road ensemble means [P05, P95]\n');
for j=1:2
    fprintf(fid,'%s, limit +/-%.0f N\n',controllerNames(j),selectedLimits(j));
    ids=summaryController==controllerNames(j);
    rows=find(ids);
    for k=1:numel(rows)
        r=rows(k);
        fprintf(fid,'  %s: %.6g [%.6g, %.6g]\n',summaryMetric(r), ...
            TRoad.Mean(r),TRoad.P05(r),TRoad.P95(r));
    end
end
fclose(fid);

disp(TDLC);
disp(TBump);
disp(TRoad);
fprintf('Reviewer 1 study completed. Outputs: %s and %s\n',dataDir,figDir);

%% Local functions
function [A,Br,BF]=intermediate_model(p)
A=[0 1 0 0 0 0;-(p.kl+p.kt)/p.mu -p.bs/p.mu p.kl/p.mu 0 0 p.bs/p.mu; ...
   0 0 0 1 0 0;p.kl/p.mi 0 -(p.kl+p.ku)/p.mi -p.ci/p.mi p.ku/p.mi p.ci/p.mi; ...
   0 0 0 0 0 1;0 p.bs/p.ms p.ku/p.ms p.ci/p.ms -p.ku/p.ms -(p.bs+p.ci)/p.ms];
Br=[0;p.kt/p.mu;0;0;0;0];
BF=[0;1/p.mu;0;0;0;-1/p.ms];
end

function [A,Br]=conventional_model(p)
A=[0 1 0 0;-(p.keq+p.kt)/p.mu -p.bs/p.mu p.keq/p.mu p.bs/p.mu; ...
   0 0 0 1;p.keq/p.ms p.bs/p.ms -p.keq/p.ms -p.bs/p.ms];
Br=[0;p.kt/p.mu;0;0];
end

function sys=passive_output_system(A,Br,p,isInt)
n=size(A,1); C=zeros(8,n); D=zeros(8,1);
if isInt
    C(1,:)=A(6,:); C(2,5)=1; C(3,1)=-1; C(3,5)=1;
    C(4,1)=p.kt; D(4)=-p.kt; C(5,3)=1;
    C(6,1)=-1; C(6,3)=1; C(7,3)=-1; C(7,5)=1;
else
    C(1,:)=A(4,:); C(2,3)=1; C(3,1)=-1; C(3,3)=1;
    C(4,1)=p.kt; D(4)=-p.kt;
end
sys=ss(A,Br,C,D);
end

function zr=bump_input(t,t0,hb,Tb)
zr=zeros(size(t));
id=t>=t0 & t<=t0+Tb;
zr(id)=hb/2*(1-cos(2*pi*(t(id)-t0)/Tb));
end

function z=iso8608_profile(t,v,seed)
rng(seed,'twister');
n=(0.01:0.005:10).'; dn=n(2)-n(1);
G=256e-6*(n/0.1).^(-2);
amp=sqrt(2*G*dn); phi=2*pi*rand(size(n));
x=(v*t).'; z=zeros(size(x));
for k=1:40:numel(n)
    id=k:min(k+39,numel(n));
    z=z+sum(amp(id).*cos(2*pi*n(id).*x+phi(id)),1);
end
z=z(:)-mean(z);
end

function [y,aux]=simulate_saturated_diagnostics(A,Br,BF,K,p,zr,t,limit)
N=numel(t); x=zeros(6,1); y=zeros(N,8); dt=t(2)-t(1);
aux.Fcmd=zeros(N,1); aux.F=zeros(N,1); aux.saturated=false(N,1);
aux.power=zeros(N,1);
for i=1:N
    Fcmd=-K*x;
    F=max(-limit,min(limit,Fcmd));
    aux.Fcmd(i)=Fcmd; aux.F(i)=F;
    aux.saturated(i)=abs(Fcmd)>limit+1e-9;
    % With BF=[0,1/mu,0,0,0,-1/ms]^T, actuator mechanical
    % power delivered to the mechanical system is F*(zu_dot-zs_dot).
    aux.power(i)=F*(x(2)-x(6));
    y(i,:)=[A(6,:)*x+BF(6)*F,x(5),x(5)-x(1), ...
        p.kt*(x(1)-zr(i)),x(3),x(3)-x(1),x(5)-x(3),F];
    if i<N
        u1=zr(i); u2=0.5*(zr(i)+zr(i+1)); u4=zr(i+1);
        ff=@(xx,u) A*xx+Br*u+BF*max(-limit,min(limit,-K*xx));
        k1=ff(x,u1); k2=ff(x+0.5*dt*k1,u2);
        k3=ff(x+0.5*dt*k2,u2); k4=ff(x+dt*k3,u4);
        x=x+dt*(k1+2*k2+2*k3+k4)/6;
    end
end
end

function m=response_metrics(y)
m.rmsAcc=rms(y(:,1)); m.peakAcc=max(abs(y(:,1)));
m.rmsTravel=rms(y(:,3)); m.maxTravel=max(abs(y(:,3)));
m.rmsTire=rms(y(:,4)); m.peakTire=max(abs(y(:,4)));
if size(y,2)>=8
    m.rmsForce=rms(y(:,8)); m.peakForce=max(abs(y(:,8)));
else
    m.rmsForce=0; m.peakForce=0;
end
end

function d=actuator_diagnostics(aux,id,dt)
sat=aux.saturated(id); P=aux.power(id);
d.saturationPercent=100*mean(sat);
d.saturationDuration=sum(sat)*dt;
d.rmsPower=rms(P);
d.peakAbsPower=max(abs(P));
end
