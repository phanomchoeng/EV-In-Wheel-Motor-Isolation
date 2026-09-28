clear; close all; clc;

rootDir = fileparts(mfilename('fullpath'));
figDir = fullfile(rootDir,'ev_iwm_figures');
dataDir = fullfile(rootDir,'ev_iwm_data');
if ~exist(figDir,'dir'), mkdir(figDir); end
if ~exist(dataDir,'dir'), mkdir(dataDir); end
set(groot,'defaultAxesFontName','Times New Roman','defaultAxesFontSize',10, ...
    'defaultAxesLineWidth',0.8,'defaultLineLineWidth',1.35,'defaultFigureColor','w');

%% Vehicle corner and road definitions
Mcorner = 400;       % kg, constant total corner mass
ms0 = 370;           % kg, centralized-drive sprung mass
mu0 = 30;            % kg, wheel-side mass without the motor
dmNom = 40;          % kg, nominal in-wheel motor and accessories
keq = 20400;         % N/m
kt = 192000;         % N/m
cs = 1400;           % N s/m

dt = 0.001; t = (0:dt:10).'; t0 = 1; hb = 0.05;
TbCases = [0.12 0.20 0.35]; speeds = [10 20 30];
n = (0.01:0.01:10).'; G = 256e-6*(n/0.1).^(-2);

pCentral = struct('ms',ms0,'mu',mu0,'kt',kt,'bs',cs,'keq',keq);
[Ac,Brc] = conventional_model(pCentral);
sysCentral = passive_output_system(Ac,Brc,pCentral,false);

pRigid = struct('ms',ms0-dmNom,'mu',mu0+dmNom,'kt',kt,'bs',cs,'keq',keq);
[Ar,Brr] = conventional_model(pRigid);
sysRigid = passive_output_system(Ar,Brr,pRigid,false);
rigidAgg = multi_scenario(sysRigid,t,t0,hb,TbCases,speeds,n,G);

%% Nominal passive search: isolate the 40 kg motor as the intermediate mass
rng(20260925,'twister'); N = 2500;
U = lhsdesign(N,2,'criterion','maximin','iterations',10);
rho = 10.^(-0.6990 + 1.3980*U(:,1));       % 0.2 to 5
ci = 10.^(2 + 1.6990*U(:,2));              % 100 to 5000 N s/m
R = nan(N,19);
for i = 1:N
    p = make_intermediate(ms0,mu0,dmNom,keq,kt,cs,rho(i),ci(i));
    [A,Br] = intermediate_model(p);
    sys = passive_output_system(A,Br,p,true);
    m = multi_scenario(sys,t,t0,hb,TbCases,speeds,n,G);
    stroke = representative_strokes(sys,t,t0,hb,0.20);
    ratios = m./rigidAgg;
    J = max(ratios);
    feasible = stroke(1) <= 0.060 && stroke(2) <= 0.060;
    if ~feasible, J = J + 10; end
    R(i,:) = [rho(i),ci(i),p.ku,p.kl,m,ratios,J,stroke];
    if mod(i,250)==0, fprintf('Nominal search %d/%d\n',i,N); end
end
R = sortrows(R,17);
best = R(1,:);
pIso = make_intermediate(ms0,mu0,dmNom,keq,kt,cs,best(1),best(2));
[Ai,Bri,BFi] = intermediate_model(pIso);
sysIso = passive_output_system(Ai,Bri,pIso,true);
isoAgg = multi_scenario(sysIso,t,t0,hb,TbCases,speeds,n,G);
writematrix(R,fullfile(dataDir,'nominal_passive_search.csv'));

%% Fixed-design and retuned behavior over representative motor masses
dmValues = (20:10:60).'; nm = numel(dmValues);
fixedRatio = nan(nm,6); tunedRatio = nan(nm,6);
tunedPar = nan(nm,5); rigidMetric = nan(nm,6);
rng(20260926,'twister'); Ntune = 800;
Ut = lhsdesign(Ntune,2,'criterion','maximin','iterations',8);
rhoT = 10.^(-0.6990 + 1.3980*Ut(:,1));
ciT = 10.^(2 + 1.6990*Ut(:,2));
for j = 1:nm
    dm = dmValues(j);
    pr = struct('ms',ms0-dm,'mu',mu0+dm,'kt',kt,'bs',cs,'keq',keq);
    [A,B] = conventional_model(pr);
    base = multi_scenario(passive_output_system(A,B,pr,false),t,t0,hb,TbCases,speeds,n,G);
    rigidMetric(j,:) = base;

    pf = make_intermediate(ms0,mu0,dm,keq,kt,cs,best(1),best(2));
    [A,B] = intermediate_model(pf);
    fixedRatio(j,:) = multi_scenario(passive_output_system(A,B,pf,true),t,t0,hb,TbCases,speeds,n,G)./base;

    bestJ = inf; bestRow = nan(1,5); bestRat = nan(1,6);
    for i = 1:Ntune
        p = make_intermediate(ms0,mu0,dm,keq,kt,cs,rhoT(i),ciT(i));
        [A,B] = intermediate_model(p); sy = passive_output_system(A,B,p,true);
        m = multi_scenario(sy,t,t0,hb,TbCases,speeds,n,G); ratios = m./base;
        stroke = representative_strokes(sy,t,t0,hb,0.20);
        J = max(ratios); if any(stroke > 0.060), J = J + 10; end
        if J < bestJ, bestJ = J; bestRow = [rhoT(i),ciT(i),p.ku,p.kl,J]; bestRat = ratios; end
    end
    tunedPar(j,:) = bestRow; tunedRatio(j,:) = bestRat;
    fprintf('Motor-mass tuning %d/%d\n',j,nm);
end

%% LQR designs for the nominal isolated-motor architecture
[Kcomfort,designComfort] = design_lqr(Ai,BFi,0.6,0.06,0.02,2000);
[Kbalanced,designBalanced] = design_lqr(Ai,BFi,1.0,0.08,0.02,1200);
sysComfort = lqr_output_system(Ai,Bri,BFi,Kcomfort,pIso);
sysBalanced = lqr_output_system(Ai,Bri,BFi,Kbalanced,pIso);

%% Figure 3: penalty caused by rigidly attaching the motor to the wheel
dmPenalty = (0:5:60).'; pen = nan(numel(dmPenalty),6);
centralAgg = multi_scenario(sysCentral,t,t0,hb,TbCases,speeds,n,G);
for j=1:numel(dmPenalty)
    dm=dmPenalty(j); p=struct('ms',ms0-dm,'mu',mu0+dm,'kt',kt,'bs',cs,'keq',keq);
    [A,B]=conventional_model(p);
    pen(j,:)=multi_scenario(passive_output_system(A,B,p,false),t,t0,hb,TbCases,speeds,n,G)./centralAgg;
end
f=figure('Position',[70 60 1160 700]); tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
a=nexttile(tl); plot(a,dmPenalty,pen(:,1:3),'-o','MarkerIndices',1:2:numel(dmPenalty)); yline(a,1,'k--'); grid(a,'on'); xlabel(a,'Motor mass transferred to the wheel (kg)'); ylabel(a,'Ratio to centralized-drive EV'); title(a,'Bump response'); legend(a,{'RMS acceleration','Peak acceleration','Maximum travel'},'Location','best');
a=nexttile(tl); plot(a,dmPenalty,pen(:,4:6),'-o','MarkerIndices',1:2:numel(dmPenalty)); yline(a,1,'k--'); grid(a,'on'); xlabel(a,'Motor mass transferred to the wheel (kg)'); ylabel(a,'Ratio to centralized-drive EV'); title(a,'Class C random road'); legend(a,{'RMS acceleration','RMS travel','RMS tire force'},'Location','best');
title(tl,'Effect of Rigid In-Wheel Motor Mass at Constant Vehicle Corner Mass','FontWeight','bold'); savefigs(f,figDir,'Figure3_Rigid_Motor_Mass_Penalty');

%% Figure 4: nominal passive design search
valid = R(:,17)<2;
f=figure('Position',[80 70 1120 470]); tl=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
a=nexttile(tl); scatter(a,R(valid,1),R(valid,17),18,log10(R(valid,2)),'filled','MarkerFaceAlpha',0.5); hold(a,'on'); scatter(a,best(1),best(17),130,'kp','filled'); grid(a,'on'); xlabel(a,'Spring stiffness ratio, \rho = k_u/k_l'); ylabel(a,'Worst normalized ratio'); title(a,'Feasible designs; color = log_{10}(c_i)'); colorbar(a);
a=nexttile(tl); scatter(a,mean(R(valid,11:13),2),mean(R(valid,14:16),2),18,R(valid,2),'filled','MarkerFaceAlpha',0.5); hold(a,'on'); scatter(a,mean(best(11:13)),mean(best(14:16)),130,'kp','filled'); grid(a,'on'); xlabel(a,'Bump index'); ylabel(a,'Random-road index'); title(a,'Tradeoff at 40 kg motor mass; color = c_i'); colorbar(a);
title(tl,'Passive Isolation Search for the Nominal In-Wheel Motor','FontWeight','bold'); savefigs(f,figDir,'Figure4_Nominal_Passive_Search');

%% Figure 5: fixed and retuned passive design over motor mass
f=figure('Position',[80 50 1180 760]); tl=tiledlayout(2,3,'TileSpacing','compact','Padding','compact'); labels={'Bump RMS acceleration','Bump peak acceleration','Maximum bump travel','Road RMS acceleration','Road RMS travel','Road RMS tire force'};
for q=1:6
    a=nexttile(tl); plot(a,dmValues,fixedRatio(:,q),'-o','DisplayName','40 kg design held fixed'); hold(a,'on'); plot(a,dmValues,tunedRatio(:,q),'-s','DisplayName','Retuned at each mass'); yline(a,1,'k--','HandleVisibility','off'); grid(a,'on'); xlabel(a,'In-wheel motor mass (kg)'); ylabel(a,'Ratio to rigid mounting'); title(a,labels{q}); if q==1, legend(a,'Location','best'); end
end
title(tl,'Passive Isolation Benefit across In-Wheel Motor Mass','FontWeight','bold'); savefigs(f,figDir,'Figure5_Motor_Mass_Robustness');

%% Figure 6: representative bump response
Tb=0.20; zrB=bump_input(t,t0,hb,Tb);
yCentral=lsim(sysCentral,zrB,t); yRigid=lsim(sysRigid,zrB,t); yIso=lsim(sysIso,zrB,t);
yBal=simulate_saturated(Ai,Bri,BFi,Kbalanced,pIso,zrB,t,1000);
f=figure('Position',[70 45 1180 760]); tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact'); idx=[1 3 4 6]; yl={'Sprung acceleration (m s^{-2})','Suspension travel (m)','Dynamic tire force (N)','Motor relative displacement z_i-z_u (m)'};
for q=1:4
    a=nexttile(tl); hold(a,'on'); grid(a,'on');
    if q<4, plot(a,t,yCentral(:,idx(q)),'k--','DisplayName','Centralized drive'); plot(a,t,yRigid(:,idx(q)),':','Color',[0.45 0.45 0.45],'DisplayName','Rigid in-wheel motor'); end
    plot(a,t,yIso(:,idx(q)),'Color',[0 0.45 0.74],'DisplayName','Passive isolated motor'); plot(a,t,yBal(:,idx(q)),'Color',[0.49 0.18 0.56],'DisplayName','Balanced LQR, +/-1000 N');
    xlim(a,[0.8 5]); xlabel(a,'Time (s)'); ylabel(a,yl{q}); if q==1, legend(a,'Location','best'); end
end
title(tl,'Response to the Representative 50 mm Bump at 40 kg Motor Mass','FontWeight','bold'); savefigs(f,figDir,'Figure6_Nominal_Bump');

%% Figure 7 and Table 4: actuator force saturation
limits=[500 1000 1500 2000 Inf]; dtR=0.002; tR=(0:dtR:60).'; zrR=iso8608_profile(tR,20,20260925); iR=tR>=5;
ctlNames=["Comfort LQR";"Balanced LQR"]; gains={Kcomfort,Kbalanced}; satRows=nan(10,8); rowCtl=strings(10,1); rowLim=nan(10,1); rr=0;
for j=1:2
 for k=1:numel(limits)
  rr=rr+1; yb=simulate_saturated(Ai,Bri,BFi,gains{j},pIso,zrB,t,limits(k)); yr=simulate_saturated(Ai,Bri,BFi,gains{j},pIso,zrR,tR,limits(k));
  mb=metrics(yb,t,t0); mr=metrics(yr(iR,:),tR(iR),5);
  rowCtl(rr)=ctlNames(j); rowLim(rr)=limits(k); satRows(rr,:)=[mb.rmsAcc mb.peakAcc 1e3*mb.maxTravel mb.peakForce mr.rmsAcc 1e3*mr.rmsTravel mr.rmsTire mr.rmsForce];
 end
end
Tsat=table(rowCtl,rowLim,satRows(:,1),satRows(:,2),satRows(:,3),satRows(:,4),satRows(:,5),satRows(:,6),satRows(:,7),satRows(:,8), ...
 'VariableNames',{'Controller','Force_limit_N','Bump_RMS_acc','Bump_peak_acc','Bump_max_travel_mm','Bump_peak_force_N','Road_RMS_acc','Road_RMS_travel_mm','Road_RMS_tire_N','Road_RMS_force_N'});
writetable(Tsat,fullfile(dataDir,'Table4_Active_Saturation.csv'));
f=figure('Position',[90 70 1120 650]); tl=tiledlayout(1,3,'TileSpacing','compact','Padding','compact'); finiteL=limits(1:4);
for q=1:3
 a=nexttile(tl); hold(a,'on'); grid(a,'on');
 for j=1:2
  ids=(j-1)*5+(1:4); if q==1, vv=satRows(ids,5); yy='Road RMS acceleration (m s^{-2})'; elseif q==2, vv=satRows(ids,6); yy='Road RMS travel (mm)'; else, vv=satRows(ids,7); yy='Road RMS tire force (N)'; end
  plot(a,finiteL,vv,'-o','DisplayName',ctlNames(j));
 end
 xlabel(a,'Actuator force limit (N)'); ylabel(a,yy); if q==1, legend(a,'Location','best'); end
end
title(tl,'Actuator Saturation for the 40 kg Isolated-Motor Architecture','FontWeight','bold'); savefigs(f,figDir,'Figure7_Active_Saturation');

%% Figure 8 and Table 5: random-road ensemble
nReal=30; ens=nan(nReal,5,3); ensNames=["Centralized drive";"Rigid in-wheel motor";"Passive isolated motor";"Comfort LQR, +/-1500 N";"Balanced LQR, +/-1000 N"];
for s=1:nReal
 zr=iso8608_profile(tR,20,20261000+s); yc=lsim(sysCentral,zr,tR); yr=lsim(sysRigid,zr,tR); yi=lsim(sysIso,zr,tR);
 ylc=simulate_saturated(Ai,Bri,BFi,Kcomfort,pIso,zr,tR,1500); ylb=simulate_saturated(Ai,Bri,BFi,Kbalanced,pIso,zr,tR,1000); YY={yc,yr,yi,ylc,ylb};
 for j=1:5, m=metrics(YY{j}(iR,:),tR(iR),5); ens(s,j,:)=[m.rmsAcc 1e3*m.rmsTravel m.rmsTire]; end
end
eRows=[]; ec=strings(15,1); em=strings(15,1); rr=0; metricNames=["RMS acceleration (m/s^2)" "RMS travel (mm)" "RMS tire force (N)"];
for j=1:5, for q=1:3, rr=rr+1; v=squeeze(ens(:,j,q)); ec(rr)=ensNames(j); em(rr)=metricNames(q); eRows(rr,:)=[mean(v) std(v) prctile(v,5) prctile(v,95)]; end, end
Tensemble=table(ec,em,eRows(:,1),eRows(:,2),eRows(:,3),eRows(:,4),'VariableNames',{'Configuration','Metric','Mean','Standard_deviation','P05','P95'});
writetable(Tensemble,fullfile(dataDir,'Table5_Random_Road_Ensemble.csv'));
f=figure('Position',[90 60 1180 720]); tl=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for q=1:3
 a=nexttile(tl); V=squeeze(ens(:,:,q)); boxchart(a,repelem((1:5).',nReal),V(:)); grid(a,'on'); set(a,'XTick',1:5,'XTickLabel',{'Central','Rigid IWM','Passive','Comfort','Balanced'}); ylabel(a,metricNames(q)); xtickangle(a,20);
end
title(tl,'Thirty Independent Class C Roads at 40 kg Motor Mass','FontWeight','bold'); savefigs(f,figDir,'Figure8_Random_Ensemble');

%% Figure 9: speed dependence
speedPlot=10:1:30; speedMetric=nan(numel(speedPlot),5,3); Sspeed={sysCentral,sysRigid,sysIso,sysComfort,sysBalanced}; speedNames={'Centralized drive','Rigid in-wheel motor','Passive isolated motor','Comfort LQR','Balanced LQR'};
for iv=1:numel(speedPlot), for j=1:5, H=squeeze(freqresp(Sspeed{j},2*pi*speedPlot(iv)*n)); speedMetric(iv,j,:)=[sqrt(trapz(n,abs(H(1,:).').^2.*G)),sqrt(trapz(n,abs(H(3,:).').^2.*G)),sqrt(trapz(n,abs(H(4,:).').^2.*G))]; end, end
f=figure('Position',[100 70 1180 620]); tl=tiledlayout(1,3,'TileSpacing','compact','Padding','compact'); yls={'RMS acceleration (m s^{-2})','RMS suspension travel (m)','RMS tire force (N)'};
for q=1:3, a=nexttile(tl); hold(a,'on'); grid(a,'on'); for j=1:5, plot(a,speedPlot,speedMetric(:,j,q),'DisplayName',speedNames{j}); end; xlabel(a,'Vehicle speed (m s^{-1})'); ylabel(a,yls{q}); title(a,yls{q}); if q==1, legend(a,'Location','best'); end, end
title(tl,'Class C Random-Road Performance for the Nominal In-Wheel Motor','FontWeight','bold'); savefigs(f,figDir,'Figure9_Speed_Sensitivity');

%% Figure 10: parameter sensitivity around the nominal isolated-motor design
nS=500; rng(20260928,'twister'); U=lhsdesign(nS,5,'criterion','maximin','iterations',10);
factor=[0.9+0.2*U(:,1),0.9+0.2*U(:,2),0.9+0.2*U(:,3),0.8+0.4*U(:,4),0.9+0.2*U(:,5)];
rob=nan(nS,3,3); stable=true(nS,3);
for i=1:nS
 p=pIso; p.ms=pIso.ms*factor(i,1); p.mi=pIso.mi*factor(i,2); p.bs=pIso.bs*factor(i,3); p.ci=pIso.ci*factor(i,4); p.kt=pIso.kt*factor(i,5);
 [A,Br,BF]=intermediate_model(p); SS={passive_output_system(A,Br,p,true),lqr_output_system(A,Br,BF,Kcomfort,p),lqr_output_system(A,Br,BF,Kbalanced,p)};
 for j=1:3
  stable(i,j)=isstable(SS{j}); H=squeeze(freqresp(SS{j},2*pi*20*n)); rob(i,j,:)=[sqrt(trapz(n,abs(H(1,:).').^2.*G)),sqrt(trapz(n,abs(H(3,:).').^2.*G)),sqrt(trapz(n,abs(H(4,:).').^2.*G))];
 end
end
sens=nan(5,3,3); for j=1:3, for q=1:3, sens(:,j,q)=corr(factor,squeeze(rob(:,j,q)),'Type','Spearman'); end, end
paramNames=["m_s";"m_i";"c_s";"c_i";"k_t"]; controllerNames=["Passive";"Comfort LQR";"Balanced LQR"];
f=figure('Position',[100 60 1180 720]); tl=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for j=1:3, a=nexttile(tl); imagesc(a,squeeze(sens(:,j,:)),[-1 1]); colorbar(a); set(a,'XTick',1:3,'XTickLabel',{'Acceleration','Travel','Tire force'},'YTick',1:5,'YTickLabel',paramNames); title(a,controllerNames(j)); for p=1:5, for q=1:3, text(a,q,p,sprintf('%.2f',sens(p,j,q)),'HorizontalAlignment','center'); end, end, end
title(tl,'Spearman Sensitivity under 500 Parameter Samples','FontWeight','bold'); savefigs(f,figDir,'Figure10_Sensitivity');

%% Tables and MAT file
Tdesign=table(["Centralized-drive EV";"Rigid in-wheel-motor EV";"Passive isolated-motor EV"], ...
 [pCentral.ms;pRigid.ms;pIso.ms],[0;0;pIso.mi],[pCentral.mu;pRigid.mu;pIso.mu], ...
 [NaN;NaN;pIso.ku],[NaN;NaN;pIso.kl],[0;0;pIso.ci],[pCentral.bs;pRigid.bs;pIso.bs], ...
 'VariableNames',{'Configuration','ms_kg','mi_kg','mu_kg','ku_Npm','kl_Npm','ci_Nspm','cs_Nspm'});
writetable(Tdesign,fullfile(dataDir,'Table1_Design_Parameters.csv'));

names=["Centralized-drive EV";"Rigid in-wheel-motor EV";"Passive isolated-motor EV"];
agg=[centralAgg;rigidAgg;isoAgg]; ratiosToRigid=agg./rigidAgg;
Tpassive=table(names,agg(:,1),agg(:,2),agg(:,3),agg(:,4),agg(:,5),agg(:,6),ratiosToRigid(:,1),ratiosToRigid(:,2),ratiosToRigid(:,3),ratiosToRigid(:,4),ratiosToRigid(:,5),ratiosToRigid(:,6), ...
 'VariableNames',{'Configuration','Bump_RMS_acc','Bump_peak_acc','Bump_max_travel','Road_RMS_acc','Road_RMS_travel','Road_RMS_tire','Ratio_bump_RMS','Ratio_bump_peak','Ratio_bump_travel','Ratio_road_acc','Ratio_road_travel','Ratio_road_tire'});
writetable(Tpassive,fullfile(dataDir,'Table2_Nominal_Passive_Performance.csv'));

Tmass=array2table([dmValues fixedRatio tunedRatio tunedPar],'VariableNames',{'Motor_mass_kg','Fixed_bump_RMS','Fixed_bump_peak','Fixed_bump_travel','Fixed_road_acc','Fixed_road_travel','Fixed_road_tire','Tuned_bump_RMS','Tuned_bump_peak','Tuned_bump_travel','Tuned_road_acc','Tuned_road_travel','Tuned_road_tire','Tuned_rho','Tuned_ci','Tuned_ku','Tuned_kl','Tuned_worst_ratio'});
writetable(Tmass,fullfile(dataDir,'Table3_Motor_Mass_Sweep.csv'));

save(fullfile(dataDir,'ev_iwm_results.mat'),'pCentral','pRigid','pIso','best','R','dmValues','fixedRatio','tunedRatio','tunedPar','Kcomfort','Kbalanced','designComfort','designBalanced','Tsat','Tensemble','ens','speedMetric','rob','stable','sens','Tdesign','Tpassive','Tmass');
fprintf('EV IWM study completed. rho=%.5f, ci=%.2f N s/m, worst nominal ratio=%.5f.\n',best(1),best(2),best(17));
fprintf('Kcomfort = [%s]\n',sprintf('%.6f ',Kcomfort));
fprintf('Kbalanced = [%s]\n',sprintf('%.6f ',Kbalanced));
fprintf('Stable samples passive/comfort/balanced: %d/%d, %d/%d, %d/%d.\n',sum(stable(:,1)),nS,sum(stable(:,2)),nS,sum(stable(:,3)),nS);

%% Local functions
function p=make_intermediate(ms0,mu0,dm,keq,kt,cs,rho,ci)
p=struct('ms',ms0-dm,'mi',dm,'mu',mu0,'kt',kt,'bs',cs,'ci',ci, ...
    'ku',keq*(1+rho),'kl',keq*(1+rho)/rho);
end
function [A,Br,BF]=intermediate_model(p)
A=[0 1 0 0 0 0;-(p.kl+p.kt)/p.mu -p.bs/p.mu p.kl/p.mu 0 0 p.bs/p.mu;0 0 0 1 0 0;p.kl/p.mi 0 -(p.kl+p.ku)/p.mi -p.ci/p.mi p.ku/p.mi p.ci/p.mi;0 0 0 0 0 1;0 p.bs/p.ms p.ku/p.ms p.ci/p.ms -p.ku/p.ms -(p.bs+p.ci)/p.ms];
Br=[0;p.kt/p.mu;0;0;0;0]; BF=[0;1/p.mu;0;0;0;-1/p.ms];
end
function [A,Br]=conventional_model(p)
A=[0 1 0 0;-(p.keq+p.kt)/p.mu -p.bs/p.mu p.keq/p.mu p.bs/p.mu;0 0 0 1;p.keq/p.ms p.bs/p.ms -p.keq/p.ms -p.bs/p.ms]; Br=[0;p.kt/p.mu;0;0];
end
function sys=passive_output_system(A,Br,p,isInt)
n=size(A,1); C=zeros(8,n); D=zeros(8,1);
if isInt
 C(1,:)=A(6,:); C(2,5)=1; C(3,1)=-1; C(3,5)=1; C(4,1)=p.kt; D(4)=-p.kt; C(5,3)=1; C(6,1)=-1; C(6,3)=1; C(7,3)=-1; C(7,5)=1;
else
 C(1,:)=A(4,:); C(2,3)=1; C(3,1)=-1; C(3,3)=1; C(4,1)=p.kt; D(4)=-p.kt;
end
sys=ss(A,Br,C,D);
end
function sys=lqr_output_system(A,Br,BF,K,p)
Acl=A-BF*K; C=zeros(8,6); D=zeros(8,1); C(1,:)=A(6,:)-BF(6)*K; C(2,5)=1; C(3,1)=-1; C(3,5)=1; C(4,1)=p.kt; D(4)=-p.kt; C(5,3)=1; C(6,1)=-1; C(6,3)=1; C(7,3)=-1; C(7,5)=1; C(8,:)=-K; sys=ss(Acl,Br,C,D);
end
function [K,d]=design_lqr(A,BF,a0,s0,t0,F0)
Ca=A(6,:); Cs=[-1 0 0 0 1 0]; Ct=[1 0 0 0 0 0]; Cz=[Ca/a0;Cs/s0;Ct/t0;zeros(1,6)]; Du=[BF(6)/a0;0;0;1/F0]; K=lqr(A,BF,Cz'*Cz,Du'*Du,Cz'*Du); d=struct('a0',a0,'s0',s0,'t0',t0,'F0',F0);
end
function zr=bump_input(t,t0,hb,Tb)
zr=zeros(size(t)); id=t>=t0&t<=t0+Tb; zr(id)=hb/2*(1-cos(2*pi*(t(id)-t0)/Tb));
end
function v=bump_triplet(sys,t,t0,hb,Tb)
y=lsim(sys,bump_input(t,t0,hb,Tb),t); i=t>=t0; v=[rms(y(i,1)) max(abs(y(i,1))) max(abs(y(i,3)))];
end
function v=road_triplet(sys,speed,n,G)
H=squeeze(freqresp(sys,2*pi*speed*n)); v=[sqrt(trapz(n,abs(H(1,:).').^2.*G)) sqrt(trapz(n,abs(H(3,:).').^2.*G)) sqrt(trapz(n,abs(H(4,:).').^2.*G))];
end
function m=multi_scenario(sys,t,t0,hb,Tbs,speeds,n,G)
b=zeros(numel(Tbs),3); for k=1:numel(Tbs), b(k,:)=bump_triplet(sys,t,t0,hb,Tbs(k)); end
r=zeros(numel(speeds),3); for k=1:numel(speeds), r(k,:)=road_triplet(sys,speeds(k),n,G); end
m=[mean(b(:,1)) mean(b(:,2)) mean(b(:,3)) mean(r(:,1)) mean(r(:,2)) mean(r(:,3))];
end
function s=representative_strokes(sys,t,t0,hb,Tb)
y=lsim(sys,bump_input(t,t0,hb,Tb),t); i=t>=t0; s=[max(abs(y(i,6))) max(abs(y(i,7)))];
end
function z=iso8608_profile(t,v,seed)
rng(seed,'twister'); n=(0.01:0.005:10).'; dn=n(2)-n(1); G=256e-6*(n/0.1).^(-2); amp=sqrt(2*G*dn); phi=2*pi*rand(size(n)); x=(v*t).'; z=zeros(size(x)); for k=1:40:numel(n), id=k:min(k+39,numel(n)); z=z+sum(amp(id).*cos(2*pi*n(id).*x+phi(id)),1); end; z=z(:)-mean(z);
end
function y=simulate_saturated(A,Br,BF,K,p,zr,t,limit)
N=numel(t); x=zeros(6,1); y=zeros(N,8); dt=t(2)-t(1);
for i=1:N
 F=max(-limit,min(limit,-K*x)); y(i,:)=[A(6,:)*x+BF(6)*F,x(5),x(5)-x(1),p.kt*(x(1)-zr(i)),x(3),x(3)-x(1),x(5)-x(3),F];
 if i<N
  u1=zr(i); u2=0.5*(zr(i)+zr(i+1)); u4=zr(i+1); ff=@(xx,u) A*xx+Br*u+BF*max(-limit,min(limit,-K*xx));
  k1=ff(x,u1); k2=ff(x+0.5*dt*k1,u2); k3=ff(x+0.5*dt*k2,u2); k4=ff(x+dt*k3,u4); x=x+dt*(k1+2*k2+2*k3+k4)/6;
 end
end
end
function m=metrics(y,t,t0)
i=t>=t0; yy=y(i,:); m.rmsAcc=rms(yy(:,1)); m.peakAcc=max(abs(yy(:,1))); m.rmsTravel=rms(yy(:,3)); m.maxTravel=max(abs(yy(:,3))); m.rmsTire=rms(yy(:,4)); m.peakTire=max(abs(yy(:,4))); if size(yy,2)>=8, m.rmsForce=rms(yy(:,8)); m.peakForce=max(abs(yy(:,8))); else, m.rmsForce=0; m.peakForce=0; end
end
function savefigs(f,d,n)
exportgraphics(f,fullfile(d,[n '.png']),'Resolution',300); exportgraphics(f,fullfile(d,[n '.pdf']),'ContentType','vector'); savefig(f,fullfile(d,[n '.fig'])); close(f);
end
