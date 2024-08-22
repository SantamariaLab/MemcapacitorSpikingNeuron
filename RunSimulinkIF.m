
%% This is the final analysis for the modeled super-capacitor
% Figure 4E
clear
open_system('fLIFv1.slx');
% %%
% % %run the simulations with R1= 1000, R2= 1e4, 2.5e4, 5e4, R3= 100. Just run the sine input tests.
% Freq=[1 1/5 1/10 1/15 1/20];%1/2 1/5 1/10 1/15 1/20]%1/2, 1/5,1/10, 1/15,1/20 1/25 1/30];
% Idc=[0.9 1.3 0.9 0.9 0.9];
% Iac=[0.02 0.1 0.1 0.1 0.1];
% %clear s
% R='[1e-2*1e3 1e3 100]';%[1000 1e3 100] 
% 
% C=['[' num2str([1*.00012 1*.0008 1*.0002]) ']'];%[.00012 .0008 .0002]
% %mdlWks = get_param('IntegrateandFire_testRC','ModelWorkspace');
% %clear(mdlWks)
% %set_param('fLIFv1/In1','LoadExternalInport','on')
% for a=1:5
%     %set_param('IntegrateandFire_testRC/AC Current Source','Frequency',num2str(Freq(a)))
%     %set_param('IntegrateandFire_testRC/Subsystem','In2',num2str(Idc(a)))
%     set_param('fLIFv1/Supercapacitor','R',R)
%     set_param('fLIFv1/Supercapacitor','C',C)
%     set_param('fLIFv1/Constant','Value',num2str(Iac(a)))
%     set_param('fLIFv1/Constant3','Value',num2str(Idc(a)))
%     %thisIdc=Idc(a);
%     %assignin('base','thisIdc',Idc(a)*[0 0; 0 0])
%     set_param('fLIFv1/Simscape Component','Frequency',num2str(Freq(a)))
%     s(a)=sim('fLIFv1.slx',240);
%     a
% end
%% s0 Run the fLIF model with 
Freq=[1 1/5 1/10 1/15 1/20];%1/2 1/5 1/10 1/15 1/20]%1/2, 1/5,1/10, 1/15,1/20 1/25 1/30];
Idc0=[0.9 0.9 0.9 0.9 0.9];
Iac0=[0.02 0.1 0.1 0.1 0.1];
%clear s
R=[1e-2*1e3 1e3 100];%[1000 1e3 100] 
C=[1*.00012 1*.0008 1*.0002];%[.00012 .0008 .0002]
backT=200;
s0=runfLIF(R,C,Freq,Idc0,Iac0,backT+11./Freq,backT);
p2a=[8 5 3 1 1];%[5 4 3 15 8 20 1 1 1];%period to analyze
[etaVF0,ampVF0,periodV0]=fracderSimscape(s0,p2a,Freq,backT);

%% S1
Freq=[1 1/5 1/10 1/15 1/20];%1/2 1/5 1/10 1/15 1/20]%1/2, 1/5,1/10, 1/15,1/20 1/25 1/30];
Idc1=[0.8 0.8 0.8 0.8 0.8];
Iac1=[0.02 0.05 0.05 0.05 0.05];
R=[5e-2*1e3 1e3 100];%[1000 1e3 100] 
C=[1*.00012 1*.0008 1*.0002];%[.00012 .0008 .0002]
backT=200;
s1=runfLIF(R,C,Freq,Idc1,Iac1,backT+11./Freq,backT);
p2a=[8 5 3 3 3 3 3 3];%[5 4 3 15 8 20 1 1 1];%period to analyze
[etaVF1,ampVF1,periodV1]=fracderSimscape(s1,p2a,Freq,backT);
%% s2
Freq=[1 1/5 1/10 1/15 1/20];%1/2 1/5 1/10 1/15 1/20]%1/2, 1/5,1/10, 1/15,1/20 1/25 1/30];
Idc2=[0.85 0.85 0.85 0.85 0.85];
Iac2=[0.02 0.05 0.05 0.05 0.05];
R=[1*1e3 1e3 100];%[1000 1e3 100] 
C=[1*.00012 1*.0008 1*.0002];%[.00012 .0008 .0002]
backT=200;
s2=runfLIF(R,C,Freq,Idc2,Iac2,backT+11./Freq,backT);
p2a=[3 3 3 3 3 3 3 3];%[5 4 3 15 8 20 1 1 1];%period to analyze
[etaVF2,ampVF2,periodV2]=fracderSimscape(s2,p2a,Freq,backT);
%% s3
Idc3=[0.85 0.85 0.85 0.85 0.85];
Iac3=[0.02 0.05 0.05 0.05 0.05];
R=[1*1e3 25*1e3 100];%[1000 1e3 100] 
C=[1*.00012 1*.0008 1*.0002];%[.00012 .0008 .0002]
backT=200;
s3=runfLIF(R,C,Freq,Idc3,Iac3,backT+11./Freq,backT);
p2a=[3 3 3 3 3 3 3 3];%[5 4 3 15 8 20 1 1 1];%period to analyze
[etaVF3,ampVF3,periodV3]=fracderSimscape(s3,p2a,Freq,backT);
%% s4
Idc4=[0.58 0.58 0.58 0.58 0.58];
Iac4=[0.01 0.02 0.05 0.05 0.05];
R=[2.5e-2*1e3 1e3 100];%[1000 1e3 100] 
C=[1*.00012 1*.0008 1*.0002];%[.00012 .0008 .0002]
backT=300;
s4=runfLIF(R,C,Freq,Idc4,Iac4,backT+11./Freq,backT);
p2a=[8 8 8 3 3];%[5 4 3 15 8 20 1 1 1];%period to analyze
[etaVF4,ampVF4,periodV4]=fracderSimscape(s4,p2a,Freq,backT);
%% s5
Idc5=[0.67 0.67 0.67 0.67 0.67];
Iac5=[0.02 0.02 0.02 0.02 0.02];
R=[3.5e-2*1e3 1e3 100];%[1000 1e3 100] 
C=[1*.00012 1*.0008 1*.0002];%[.00012 .0008 .0002]
backT=250;
s5=runfLIF(R,C,Freq,Idc5,Iac5,backT+11./Freq,backT);
p2a=[8 8 8 8 8];%[5 4 3 15 8 20 1 1 1];%period to analyze
[etaVF5,ampVF5,periodV5]=fracderSimscape(s5,p2a,Freq,backT);
%%
figure(1)
clf
h=plot(1./Freq,etaVF0*2/pi,1./Freq,etaVF1*2/pi,1./Freq,etaVF2*2/pi,...
    1./Freq,etaVF3*2/pi,1./Freq,etaVF4*2/pi,1./Freq,etaVF5*2/pi);
ylabel('Phase (rad)')
xlabel('Period (sec)')
h(1).Marker='sq';h(1).MarkerSize=10;h(1).Color=[1 0 0];h(1).MarkerFaceColor=[1 0 0];
h(2).Marker='^';h(2).MarkerSize=10;h(2).Color=[0.4 0 0];h(2).MarkerFaceColor=[0.4 0 0];
h(3).Marker='v';h(3).MarkerSize=10;h(3).Color=[0.2 0 0];h(3).MarkerFaceColor=[0.2 0 0];
h(4).Marker='>';h(4).MarkerSize=10;h(4).Color=[0.1 0 0];h(4).MarkerFaceColor=[0.1 0 0];
h(5).Marker='<';h(5).MarkerSize=10;h(5).Color=[0.8 0 0];h(5).MarkerFaceColor=[0.8 0 0];
h(6).Marker='hexagram';h(6).MarkerSize=10;h(6).Color=[0.6 0 0];h(6).MarkerFaceColor=[0.6 0 0];

ylim([0 0.6])
box  off
h2=h.Parent;
h2.TickLength=[0.05 0.025];
h2.FontSize=10;
h2.YLabel.FontSize=12;
h2.XLabel.FontSize=12;
figh=gcf;
set(0, 'DefaultFigureRenderer', 'painters');
set(figh,'WindowStyle','normal')
set(figh, 'PaperUnits', 'inches','Renderer','painters');
set(figh,'Units','inches','Position',[1 1 2 2])
exportgraphics(figh,'spicecircuitfLIF_etafLIF.eps')
set(figh,'WindowStyle','docked')



figure(2)
clf
g0=ampVF0./Iac0;
g1=ampVF1./Iac1;
g2=ampVF2./Iac2;
g3=ampVF3./Iac3;
g4=ampVF4./Iac4;
g5=ampVF5./Iac5;
gN0=g0./g0(1);
gN1=g1./g1(1);
gN2=g2./g2(1);
gN3=g3./g3(1);
gN4=g4./g4(1);
gN5=g5./g5(1);
h=loglog(1./Freq,gN0,...
    1./Freq,gN1,...
    1./Freq,gN2,...
    1./Freq,gN3,...
    1./Freq,gN4,...
    1./Freq,gN5,'color',[0 0 0]);
h(1).Marker='sq';h(1).MarkerSize=10;h(1).Color=[1 0 0];h(1).MarkerFaceColor=[1 0 0];
h(2).Marker='^';h(2).MarkerSize=10;h(2).Color=[0.4 0 0];h(2).MarkerFaceColor=[.4 0 0];
h(3).Marker='v';h(3).MarkerSize=10;h(3).Color=[0.2 0 0];h(3).MarkerFaceColor=[.2 0 0];
h(4).Marker='>';h(4).MarkerSize=10;h(4).Color=[0.1 0 0];h(4).MarkerFaceColor=[.1 0 0];
h(5).Marker='<';h(5).MarkerSize=10;h(5).Color=[0.8 0 0];h(5).MarkerFaceColor=[0.8 0 0];
h(6).Marker='hexagram';h(6).MarkerSize=10;h(6).Color=[0.6 0 0];h(6).MarkerFaceColor=[0.6 0 0];
ylabel('Norm. Gain')
xlabel('Period (sec)')
box  off
h2=h.Parent;
h2.TickLength=[0.05 0.025];
h2.FontSize=10;
h2.YLabel.FontSize=12;
h2.XLabel.FontSize=12;

figh=gcf;
set(0, 'DefaultFigureRenderer', 'painters');
set(figh,'WindowStyle','normal')
set(figh, 'PaperUnits', 'inches','Renderer','painters');
set(figh,'Units','inches','Position',[1 1 2 2])
exportgraphics(figh,'spicecircuitfLIF_Gain.eps')
set(figh,'WindowStyle','docked')
% gf0=fit(1./Freq',gN0','poly1')
% gf1=fit(1./Freq',gN1','poly1')
% gf2=fit(1./Freq',gN2','poly1')
% gf3=fit(1./Freq',gN3','poly1')
% gf4=fit(1./Freq',gN4','poly1')
% gf5=fit(1./Freq',gN5','poly1')

%% compare to super-capacitor props
% this will make Figure 4E bottom
% Figure 4A and 4B
open_system('RC')
%%
RM(1,:)=[1e-2*1e3 1e3 100];%[1000 1e3 100] 
RM(2,:)=[5e-2*1e3 1e3 100];%[1000 1e3 100] 
RM(3,:)=[1*1e3 1e3 100];%[1000 1e3 100] 
RM(4,:)=[1*1e3 25*1e3 100];%[1000 1e3 100] 
RM(5,:)=[2.5e-2*1e3 1e3 100];%[1000 1e3 100] 
RM(6,:)=[3.5e-2*1e3 1e3 100];%[1000 1e3 100]
C=[1*.00012 1*.0008 1*.0002];%[.00012 .0008 .0002]
x0=[27 1 0.2];
Iccf=0.05*1e-3;
options = optimoptions('lsqnonlin');
options.Algorithm = 'levenberg-marquardt';
options.TolFun=1.5e-100; %Magnitude of search direction was smaller than the specified tolerance.
options.TolX=1e-100;  % Change in x was less than the specified tolerance.
options.MaxIter = 100; %Number of iterations exceeded options.MaxIter
options.MaxFunEvals=1000;
clear x
figure(3)
clf
for a=1:6
    set_param('RC/Supercapacitor','R',['[' num2str(RM(a,:)) ']']);
    rcs(a)=sim('RC.slx',15);
    thisD=rcs(a).simout1.Data;
    thisT=rcs(a).simout1.Time;
    r2fit=find((thisT>=10).*(thisT<13));
    nt2f=thisT([r2fit]);
    nt2f=nt2f-nt2f(1);
    ws2f=thisD([r2fit]);
    ws2f=ws2f-ws2f(1);
    Iccf=str2double(get_param('RC/Constant3','Value'))*1.e-3;
    x0=[0.8 0.2];
    fun = @(x)x(2).*((nt2f.^x(1))/(gamma(1+x(1))))-ws2f;
    [x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(fun,x0,[],[],options);
    ci=nlparci(x,residual,'Jacobian',jacobian);
    ci2(a,:)=ci(1,:); %residuals
    etaRC(a)=x(1); %frac exp
    nt2=nt2f;%0:dt:time2(end)+0;
    v2=x(2).*((nt2.^x(1))/(gamma(1+x(1))));
    loglog(nt2f,ws2f,'LineWidth',2)
    hold on
    loglog(nt2,v2,'k--','LineWidth',2)
    %xlim([1 500])
end

% the average fractional order of the spiking circuit
etaVFM=[mean(etaVF0) mean(etaVF1) mean(etaVF2) mean(etaVF3) mean(etaVF4) mean(etaVF5)];
figure(5)
clf
f3=fit(etaRC([1 2 3 4 5 6])',etaVFM([1 2 3 4 5 6])','poly1');
%h=plot(etaRC([1 2 3 4 5 6]),etaVFM([1 2 3 4 5 6]),'o','markersize',10,'color',[0 0 0]);
h=plot(etaRC([1 ]),etaVFM([1 ]),'sq','markersize',10,'color',[1 0 0],'markerfacecolor',[1 0 0]);
hold on
h=plot(etaRC([2 ]),etaVFM([2 ]),'^','markersize',10,'color',[0.4 0 0],'markerfacecolor',[0.4 0 0]);
h=plot(etaRC([3]),etaVFM([3]),'v','markersize',10,'color',[0.2 0 0],'markerfacecolor',[0.2 0 0]);
h=plot(etaRC([4]),etaVFM([4]),'>','markersize',10,'color',[0.1 0 0],'markerfacecolor',[0.1 0 0]);
h=plot(etaRC([5]),etaVFM([5]),'<','markersize',10,'color',[0.8 0 0],'markerfacecolor',[0.8 0 0]);
h=plot(etaRC([6]),etaVFM([6]),'hexagram','markersize',10,'color',[0.6 0 0],'markerfacecolor',[0.6 0 0]);
xlabel('Frac. Ord. Super Cap.')
ylabel('Frac. Ord. fLIF')

hold on
h=plot([0:1],f3(0:1),'k')
xlim([0.35 0.55])
ylim([0 1])
box off
h2=h.Parent;
h2.TickLength=[0.05 0.025];
h2.FontSize=10;
h2.YLabel.FontSize=12;
h2.XLabel.FontSize=12;

figh=gcf;
set(0, 'DefaultFigureRenderer', 'painters');
set(figh,'WindowStyle','normal')
set(figh, 'PaperUnits', 'inches','Renderer','painters');
set(figh,'Units','inches','Position',[1 1 2 2])
exportgraphics(figh,'spicecircuitfLIF_RCetaVsfLIFeta.eps')
set(figh,'WindowStyle','docked')

%% run sweeps for super-capacitor orders
% Figure 4B
R1_o=1e3;
r1=[0:10:100 200:100:1000 2000:1000:10000];
for a=1:length(r1)
    Rin=[r1(a) 1e3 100];
    [sr1,etar1V(a)]=runSCSuperCap(Rin,15);
end
r2=[0:10:100 200:100:1000 2000:1000:10000];
for a=1:length(r1)
    Rin=[1e3 r2(a) 100];
    [sr2,etar2V(a)]=runSCSuperCap(Rin,15);
end
r3=[0:10:100 200:100:1000 2000:1000:10000];
for a=1:length(r1)
    Rin=[1e3 1e3 r3(a)];
    [sr3,etar3V(a)]=runSCSuperCap(Rin,15);
end


clf
h=plot(r1,etar1V,r2,etar2V,r3,etar3V);
h(1).Marker='o';h(1).MarkerSize=5;
h(2).Marker='o';h(2).MarkerSize=5;
h(3).Marker='o';h(3).MarkerSize=5;
xlabel('Resistance (Ohms)')
ylabel('Frac. Ord. Super Cap.')
xlim([0 5000])
box off
h2=h.Parent;
h2.TickLength=[0.05 0.025];
h2.FontSize=10;
h2.YLabel.FontSize=12;
h2.XLabel.FontSize=12;
figh=gcf;
set(0, 'DefaultFigureRenderer', 'painters');
set(figh,'WindowStyle','normal')
set(figh, 'PaperUnits', 'inches','Renderer','painters');
set(figh,'Units','inches','Position',[1 1 2 2])
exportgraphics(figh,'spicecircuitfLIF_RCetaVsR1R2R3.eps')
set(figh,'WindowStyle','docked')
%% make the sample figures of sine input
t=0:.01:10;
x=sin(2*pi*t/1);
clf
h=plot(t,x);
box off
axis off
figh=gcf;
set(0, 'DefaultFigureRenderer', 'painters');
set(figh,'WindowStyle','normal')
set(figh, 'PaperUnits', 'inches','Renderer','painters');
set(figh,'Units','inches','Position',[1 1 2 0.5])
exportgraphics(figh,'spicecircuitfLIF_sinesample.eps')
set(figh,'WindowStyle','docked')

t=0:.01:10;
x=sign(sin(2*pi*t/1));
clf
h=plot(t,x,'k');
box off
axis off
figh=gcf;
set(0, 'DefaultFigureRenderer', 'painters');
set(figh,'WindowStyle','normal')
set(figh, 'PaperUnits', 'inches','Renderer','painters');
set(figh,'Units','inches','Position',[1 1 2 0.5])
exportgraphics(figh,'spicecircuitfLIF_squaresample.eps')
set(figh,'WindowStyle','docked')

