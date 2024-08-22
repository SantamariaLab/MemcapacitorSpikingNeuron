
%%
%Analysis of the fractional and integratl LIF circuit
%%

clear 
load dIF
%%
%dIF % has firing rate in response to long lasting constant input current.
% dIFint spiking vcontant input current in the original classical circuit
%dIFSq response to square input
%[ISI,IFR,time] =ISIIFRICircuit(dIF.Exp(1));


%find a nice spike
thisV=-dIF.Exp(6).Data.Vout(2e5:3e5);
thisSp=find(thisV>3.4);
dt=1e-4;
win2p=[-1e-3/dt:10e-3/dt];
thist=1000*win2p*dt;
theSp=thisV(thisSp(1)+win2p);

%load the integer LIF
load dIFint
thisVi=-dIFint.Exp(6).Data.Vout;
thisSpi=find(thisVi>2.4);
dt=1e-4;
thist=1000*win2p*dt;
theSpi=thisVi(thisSpi(1)+win2p);

clf
plot(thist,theSp,'color',[0 0 0])
hold on
plot(thist,theSpi,'color',[0.5 0.5 0.5])
box off
axis square
ylim([1 3.7])
h=gca;
h2=h.Parent;
set(h2,'PaperUnits','inches')
set(h2,'WindowStyle','normal');
set(h2,'Units','inch')
set(h2,'Position',[4 4 2 2]);
exportgraphics(h,'SpikeExample_fLIF.eps')


%% IF plot for both of them
% dIF CV is wrong. It should be 1e-3 (1 mA per V)
% dIFint it should be 1e-3
clf
for winS=[0.5 1 10]
    %winS=0.5;%sec
    fr=0;fri=0;isi=0;
    iV=0;%[0:0.25:2];
    vV=0;
    for a=1:8
        thisV=-dIF.Exp(a).Data.Vout(:);
        iV=[iV dIF.Exp(a).Data.Vin(1)*1e-3]; %Vin is the voltage command 1 mA/V in the current source
        vV=[vV dIF.Exp(a).Data.Vin(1)];
        thisSp=find(diff(thisV)<-0.5);
        frW=thisSp(1)+[0 winS/dt];
        sp2a=(thisSp>=frW(1)).*(thisSp<=frW(2));
        fr=[fr sum(sp2a)./winS];
        isi=[isi std(diff(sp2a))./sqrt(sum(sp2a))];
        
        
        thisVi=-dIFint.Exp(a).Data.Vout(:);
        thisSpi=find(diff(thisVi)<-0.1);
        frWi=thisSpi(1)+[0 winS/dt];
        fri=[fri sum((thisSpi>=frWi(1)).*(thisSpi<=frWi(2)))./winS];
        figure(1)
        plot(thisV(1:end))
        drawnow
        input([ num2str(vV(a+1)) ' r'])
    end
    figure(2)
    plot(iV*1e3, fr,'color',[0 0 0],'LineWidth',2);
    hold on
    plot(iV*1e3,fri,'color',[0.5 0.5 0.5],'LineWidth',2)
end
box off
axis square
h=gca;
h2=h.Parent;
set(h2,'PaperUnits','inches')
set(h2,'WindowStyle','normal');
set(h2,'Units','inch')
set(h2,'Position',[4 4 2 2]);
exportgraphics(h,'FvsI_LIFandfLIF.eps')
set(h2,'WindowStyle','docked')

%% ISI 
clf
for winS=[10]
    %winS=0.5;%sec
    fr=[];fri=[];isi=[];
    iV=[0:0.25:2];
    iV=[];
    for a=7:8
        thisV=-dIF.Exp(a).Data.Vout(:);
        thisSp=find(diff(thisV)<-0.4);
        iV=[iV dIF.Exp(a).Data.Vin(1)*1e-3]; 
        frW=thisSp(1)+[0 winS/dt];
        sp2a=(thisSp>=frW(1)).*(thisSp<=frW(2));
        selSp=thisSp(logical(sp2a));
        isi{a}=[ diff(thisSp(logical(sp2a)))*dt*1000];
        

        thisVi=-dIFint.Exp(a).Data.Vout(:);
        thisSpi=find(diff(thisVi)<-0.1);
        frWi=thisSpi(1)+[0 winS/dt];
        sp2ai=(thisSpi>=frWi(1)).*(thisSpi<=frWi(2));
        isii{a}= diff(sp2ai)*dt*1000;
        
        loglog(selSp(1:5:end-1)*dt-selSp(1)*dt,isi{a}(1:5:end),'color',[0 0 0])
%         [N,edges]=histcounts(isi{a},30);
%         loglog(edges(1:end-1),N)
         hold on
       % input('r')
    end  
    
end
box off
axis square
h=gca;
h2=h.Parent;
set(h2,'PaperUnits','inches')
set(h2,'WindowStyle','normal');
set(h2,'Units','inch')
set(h2,'Position',[4 4 2 2]);
exportgraphics(h,'loglogISI_examples.eps')



%% examples of time to first spike
% the current source was set at 0.1 mA/V for Vin
a=8;
thisV=-dIF.Exp(a).Data.Vout(:);
thisSp=find(diff(thisV)<-0.5);
frW=thisSp(1);


thisVi=-dIFint.Exp(a).Data.Vout(:);
thisSpi=find(diff(thisVi)<-0.1);
clf
plot([1:thisSp(30)]*dt,thisV(1:thisSp(30)),'color',[0 0 0],'LineWidth',1)
hold on
%plot(thisVi(1:thisSpi(3)),'color',[0.5 0.5 0.5],'LineWidth',2);
ylim([0 2])
xlim([2.5 3.5])
box off
h=gca;
h2=h.Parent;
set(h2,'WindowStyle','normal');
set(h2,'Units','inch')
set(h2,'Position',[4 4 2 1]);
exportgraphics(h,'TimetoSpike.eps')

%% analyze square input
% Data.CV is the convertion factor in the current supplu so I=CV*Vin

clear
load dIFSq
out=dIFSq;
clear dIFSq
%%
p2a=1*[3 5 9 15 20 27  4];%period to analyze
spdth=[0.3 0.2 0.3 0.5 0.5 0.3 0.3]
fws=[1000 1000 1000 500 400 300 1000];
trimW=[850 600 600 300 200 280 1500];
trimWU=[100 110 350 500 350 400 130];
trimE=[600 600 600 600 600 600 600]
   
f=fittype('a*(1-exp(-x./b))');
g=fitoptions(f);
%=0;%SimProp.dt;
%tcL=[];tcU=[];
clf
c=1;
for a=1:7
    dt=out.Exp(a).Data.dt;
    %thisInp=0.01*-square(2*pi*out.Exp(a).Data.t/wl)'+0.15;
    thisVL=out.Exp(a).Data.Vout(0/dt+1:end);
    %thisS=SimProp.IinjDC(preT(a)/SimProp.dt+1:end);
    thist=out.Exp(a).Data.t;%0:SimProp.dt:numel(thisVL);
    spP=1*find(diff(thisVL>1.95)==-1);%find(diff(thisVL>-40)==-1);
    sptA=spP*dt;
    ifrA=1./diff(sptA);
    wlV=1./(out.Exp(a).Data.Frequency);
    fri=zeros(size(thist));
    %fri(spP)=ifrA;
    for b=1:length(ifrA)-1
        fri(spP(b):spP(b+1))=ifrA(b);
    end
    fri=smoothdata(fri,'movmean',fws(a));
    %extract 0.5 cycle
    clf
    w2e=(thist>=(wlV*p2a(a)+trimW(a)*dt)).*(thist<(wlV*(p2a(a)+0.5)-trimE(end)*dt));
    thisFr=fri(logical(w2e)); 
    thisFrR=(thisFr-thisFr(1));
    thistp=(0:length(thisFrR)-1).*dt;
    %thistp=thistp-thistp(1);
    g.StartPoint=[max(thisFrR) wlV/2 ];
    f0=fit(thistp(logical((thistp>0).*(thistp<16)))',thisFrR(logical((thistp>0).*(thistp<16))),f,g);
    plot(thistp,thisFrR,thistp,f0(thistp))
    tL(a)=f0.b;
    hold
    
    w2e=(thist>=(wlV*(p2a(a)+0.5+trimWU(a)*dt))).*(thist<(wlV*(p2a(a)+1.0)-trimE(end)*dt));
    thisFr=fri(logical(w2e)); 
    thisFrR=-(thisFr-thisFr(1));
    thistp=(0:length(thisFrR)-1).*dt;
    g.StartPoint=[max(thisFrR) wlV/2 ];
    f0=fit(thistp(logical(thistp<16))',thisFrR(logical(thistp<16)),f,g);
    plot(thistp,thisFrR,thistp,f0(thistp))
    tU(a)=f0.b;
    periodV(a)=wlV;
    %input('r')
    c=c+1;
end
%% Figure 3B
clf
[~,xp]=sort(periodV);

h=plot(periodV(xp),tU(xp),'-k.','linewidth',2,'MarkerSize',20);
hold on
h=plot(periodV(xp),tL(xp),'--k.','LineWidth',2,'MarkerSize',20);
h2=h.Parent;
h2.TickLength=[0.05 0.025];
box off
ylim([0 2])
% load the simulated data 
sqfLIFmodel=load('model_fLIF_SquaredInput.mat');
h=plot(sqfLIFmodel.wlV(1:7)/1000,sqfLIFmodel.tcL,'-','linewidth',2,'Color',[0.8 0.8 0.8]);
h=plot(sqfLIFmodel.wlV(1:7)/1000,sqfLIFmodel.tcU,'--k','linewidth',2,'Color',[0.8 0.8 0.8]);
h2=h.Parent;
h2.TickLength=[0.05 0.25];
box off

%% export Fig 3B

figh=gcf;
set(0, 'DefaultFigureRenderer', 'painters');
set(figh,'WindowStyle','normal')
set(figh, 'PaperUnits', 'inches','Renderer','painters');
set(figh,'Units','inches','Position',[1 1 2 2])
exportgraphics(figh,'circuitfLIF_SquaredInput_circuitVsimulations.eps')
set(figh,'WindowStyle','docked')

%% Figure 3A example square input
clf
a=1
dt=out.Exp(a).Data.dt;

thist=out.Exp(a).Data.t;
    wl=1/out.Exp(a).Data.Frequency;
    thisVL=out.Exp(a).Data.Vout(1:end);
    thisS=out.Exp(a).Data.Vin(1:end)*out.Exp(a).Data.CV;
    spP=1*find(diff(thisVL>1.80)==-1);
    sptA=spP*dt;
    ifrA=1./diff(sptA);
    fri=zeros(size(thist));
    %fri(spP)=ifrA;
    for b=1:length(ifrA)-1
        fri(spP(b):spP(b+1))=ifrA(b);
    end
clf
%tiledlayout(2,1)
%nexttile
colororder([0 0 0; 0.5 0.5 0.5])
yyaxis left
plot(thist-15,smoothdata(fri,'movmean',1000),'k','LineWidth',2)
xlim([0 20])

yyaxis right
plot(thist-15,thisS,'color',[0.5 0.5 0.5],'LineWidth',2)
xlim([0 20])
%ylim([1 2.5])
box off
figh=gcf;

set(0, 'DefaultFigureRenderer', 'painters');
set(figh,'WindowStyle','normal')
set(figh, 'PaperUnits', 'inches','Renderer','painters');
set(figh,'Units','inches','Position',[1 1 2 2])
exportgraphics(figh,'fLIFcircuit_SquareInput_Example.eps')
set(figh,'WindowStyle','docked')



%% analyze sine input
clear 
load dIFSineFinal
thisTest=dIFSineFinal;
%%
p2a=1*[40 30 1 4 15 20  20];%period to analyze
f=fittype('a*sin(2*pi*x/l+eta*pi/2)','problem','l');
g=fitoptions(f);
preT=[0 0 0 0 0 0 0 ];%5000
dt=thisTest(1).Exp(1).Data.dt;

%=0;%SimProp.dt;
%tcL=[];tcU=[];
clf
c=1;
st2a=1;
for a=1:8%1:length(wlV)
    wl=1/thisTest(st2a).Exp(a).Data.Frequency;
    p2a(a)=round(30/wl);

    thisV=(-thisTest(st2a).Exp(a).Data.Vout);
    thisVin=thisTest(st2a).Exp(a).Data.Vin;
    thisVin=thisVin-mean(thisVin);
    thist=thisTest(st2a).Exp(a).Data.t;
    thisSp=find(diff((diff(thisV)<-0.5))==1);
    thisSp=find(diff((diff(thisV)<-0.5))==1);
    isi=diff(thisSp)*dt;
    fr=1./isi;
    dfr=diff(fr);
    ifr2=zeros(size(thist));
    ifr2(thisSp(2:end-1))=dfr;
    ifr2(thisSp(1))=fr(1);
    ifr=cumsum(ifr2);

    w2e=(thist>=(wl*p2a(a)-dt)).*(thist<(wl*(p2a(a)+1.0)-dt));
    thisVin=thisVin(logical(w2e));
    thisFr=ifr(logical(w2e)); 
    thisFrR=thisFr-mean(thisFr);
    thist2=(0:length(thisFrR)-1).*dt;
    g.StartPoint=[(max(thisFrR)-min(thisFrR))/2 0.4];
    f0=fit(thist2(1:end-1)',thisFrR(1:end-1),f,g,'problem',wl);
    plot(2*pi*thist2/max(thist2),thisFrR,2*pi*thist2./max(thist2),f0(thist2),...
        2*pi*thist2./max(thist2),max(thisFrR)*thisVin/max(thisVin))
    drawnow
    wlV(a)=wl;
    ampVF(a)=f0.a;
    etaVF(a)=f0.eta;
    %hold on
    c=c+1
    %input('r')
end
%% Figure 3D
ampV=0.5;
clf
%open fLIFmodel_gain.fig
figure(2)
%plot the gain
colororder([0 0 1; 00 .5 0])
yyaxis left
gainV=ampVF./ampV;
h=loglog(wlV(1:8),gainV./max(gainV),'marker','.','markersize',20,'color',[0 0 0.5],'LineWidth',1)
hold on
fg=fit(log(wlV)',log(gainV'./max(gainV)),'poly1');
loglog(exp(log(wlV)),exp(fg(log(wlV))),'k--');
h2=h.Parent;
h2.TickLength=[0.05 0.025]
box off
ylim([.25 1])
xlim([0.5 40])
%ylabel('Norm. Gain')

%plot the phase lag
yyaxis right
h=semilogx(wlV(1:8),etaVF*2/pi,'marker','.','markersize',20,'color',[0  0.5 0],'LineWidth',1);
h2=h.Parent;
h2.TickLength=[0.05 0.025];
box off
ylim([0 1])
%ylabel('Phase (Gradians)')


figh=gcf;
set(0, 'DefaultFigureRenderer', 'painters');
set(figh,'WindowStyle','normal')
set(figh, 'PaperUnits', 'inches','Renderer','painters');
set(figh,'Units','inches','Position',[1 1 2 2])
exportgraphics(figh,'fLIFcircuit_FracDiff_gain.eps')
set(figh,'WindowStyle','docked')



% %figure (3)
% h=semilogx(wlV(1:8),etaVF,'marker','.','markersize',20,'color',[0 0 0],'LineWidth',2)
% h2=h.Parent;
% h2.TickLength=[0.05 0.025];
% box off
% ylim([0 1])
% 
% figh=gcf;
% set(0, 'DefaultFigureRenderer', 'painters');
% set(figh,'WindowStyle','normal')
% set(figh, 'PaperUnits', 'inches','Renderer','painters');
% set(figh,'Units','inches','Position',[1 1 2 2])
% exportgraphics(figh,'fLIFcircuit_FracDiff_phase.eps')
% set(figh,'WindowStyle','docked')

%% make an example of the sine input
% Data.CV is the convertion factor in the current supplu so I=CV*Vin
clf
a=3
st2a=1;
wl=1/thisTest(st2a).Exp(a).Data.Frequency;
p2a(a)=round(30/wl);
dt=thisTest(1).Exp(1).Data.dt;
thisV=(-thisTest(st2a).Exp(a).Data.Vout);
thisVin=thisTest(st2a).Exp(a).Data.Vin;
thisVin=thisVin-mean(thisVin);
thist=thisTest(st2a).Exp(a).Data.t;
thisSp=find(diff((diff(thisV)<-0.5))==1);
thisSp=find(diff((diff(thisV)<-0.5))==1);
isi=diff(thisSp)*dt;
fr=1./isi;
dfr=diff(fr);
ifr2=zeros(size(thist));
ifr2(thisSp(2:end-1))=dfr;
ifr2(thisSp(1))=fr(1);
ifr=cumsum(ifr2);

w2e=(thist>=(wl*p2a(a)-dt)).*(thist<(wl*(p2a(a)+1.0)-dt));
thisVin=thisVin(logical(w2e));
thisFr=ifr(logical(w2e));
thisFrR=thisFr-mean(thisFr);
thist2=(0:length(thisFrR)-1).*dt;
g.StartPoint=[(max(thisFrR)-min(thisFrR))/2 0.4];
f0=fit(thist2(1:end-1)',thisFrR(1:end-1),f,g,'problem',wl);

colororder([0 0 0; 0.5 0.5 0.5])
yyaxis left
h=plot(2*pi*thist2/max(thist2),thisFrR,2*pi*thist2./max(thist2),f0(thist2));
yyaxis right
h=plot(2*pi*thist2/max(thist2),thisVin,'color',[0.5 0.5 0.5],'LineWidth',2)
h2=h.Parent;
h2.XTick=[0 pi 2*pi]
h2.XTickLabel={'0','\pi','2\pi'}
xlim([0 2*pi]);
box off
figh=gcf;

set(0, 'DefaultFigureRenderer', 'painters');
set(figh,'WindowStyle','normal')
set(figh, 'PaperUnits', 'inches','Renderer','painters');
set(figh,'Units','inches','Position',[1 1 2 2])
exportgraphics(figh,'fLIFcircuit_FracDiff_Example.eps')
set(figh,'WindowStyle','docked')



%% Analyze the pink noise - Figure 4
clear
load dpnIF2.mat

%%
addpath('C:\Users\jfu936\OneDrive - University of Texas at San Antonio\Manuscripts\2020\FracCircuit\Hodgkin and Huxley\Data and code\Figure 1  ')
dt=dpnIF2(1).Exp(1).Data.dt;
a=1
clf
clear avPN
slopeV=[];
whiteIndex=[];
begf=1;
endf=8;
areaT=0.9^1; 
xgauss=-0.1:dt:0.1;
    sgauss=1*0.002;%I've used 0.002 for most analysis
    mygauss=dt/sqrt(2*pi*sgauss.^2).*exp(-0.5*xgauss.^2./sgauss.^2);%
for a=1:5
    %thisFTout=myLIFpowerplot(dpnIF2(1).Exp(a),'all');
    thisV=dpnIF2(1).Exp(a).Data.Vout;
    thist=dpnIF2(1).Exp(a).Data.t;
    thisSp=find(diff((diff(thisV)<-0.5))==1);
    isi=diff(thisSp)*dt;
    fr=1./isi;
    dfr=diff(fr);
    ifr2=zeros(size(thist));
    ifr2(thisSp(2:end-1))=dfr;
    ifr2(thisSp(1))=fr(1);
    ifr=cumsum(ifr2);
    dummy1=conv(ifr,mygauss,'same');%conv ifr to smooth edges avoid artifacts
    %thisFTout=getFFT_LIF(dummy1,1/dt,10);
    c=1;binV=[];fitV=[];outM=[];
    for c=1:8
        fr2a=dummy1(c*10/dt+1+[0:10/dt]);
        out=getFFT_LIF(fr2a-mean(fr2a),1./dt,10);
        binV(:,c)=out.BinData(:,2);
        thisFitout=fit(log10(out.BinData(11:28,1)),...
            log10(out.BinData(11:28,2)./out.BinData(28,2)),'a.*(x)+b');
        thisFitouM{c}=thisFitout;
        fitV(c)=thisFitout.a;
    end
    binVm=mean(binV,2);
    %why did I comment this next line?
    thisFTout.BinData=[out.BinData(:,1) binVm];
    mfitV=mean(fitV);
    sefitV=std(fitV)./sqrt(length(fitV));
    sertVM(a)=sefitV;

    %fitting to the average trace
    % thisFitout=fit(log10(out.BinData(10:17,1)),...
    %     log10(binVm(10:17,1)./binVm(10,1)),'a.*(x)+b');
    % mfitV=thisFitout.a;

    %analyze the input signal
    thisVin1=dpnIF2(1).Exp(a).Data.Vin;
    thisFTin=getFFT_LIF(thisVin1,1/dt,10);
    thisFitin=fit(log10(thisFTin.BinData(10:28,1)),log10(thisFTin.BinData(10:28,2)./thisFTin.BinData(10,2)),'a*(x)+b');
    %integrate the area under the curve to calculate the white index. The
    %area is 1 log unit in Hz ([0.1 1]* 1 log unit in power [0.1 1]
    outInt=trapz(log10(thisFTout.BinData(10:18,1)), log10(binVm(10:18,1)./binVm(10,1)));
    inInt=trapz(log10(thisFTin.BinData(10:18,1)),log10(thisFTin.BinData(10:18,2)./thisFTin.BinData(10,2)));
    whiteIndex=[whiteIndex; 1+outInt/(1*1) 1+inInt];
    %slopeV=[slopeV; thisFitin.a thisFitout.a]
    %inputfit outputfit
    slopeV=[slopeV; thisFitin.a mfitV]
    

    clf
    loglog(thisFTout.BinData(10:29,1),binVm(10:29,1)./binVm(10,1),'o-k')
    %loglog(thisFTout.BinData(10:29,1),thisFTout.BinData(10:29,2)./thisFTout.BinData(10,2),'o-k')
    hold on
    %the mean fit. 
    loglog(0.1:0.1:1,0.1.^-mfitV.*([0.1:0.1:1]).^mfitV,'--k')

    %loglog(thisFTin.BinData(10:28,1),thisFTin.BinData(10:28,2)./thisFTin.BinData(10,2),'.-k')
    %hold on
    %the fit to the input signal. 
    loglog(thisFTin.BinData(10:19,1),10.^(thisFitin.b).*(thisFTin.BinData(10:19,1)).^thisFitin.a,'k')

    xlim([0.1 1])
    ylim([0.1 3])
    suffix={'002','004','006','008','100'}
    if sum(a==[1 2 3 ])
        xlabel('Frequency (Hz)')
        ylabel('Power (norm.)')
        set(gca,'clipping','on')
        box off
        h=gca;
        h.FontSize=10;
        h.YLabel.FontSize=12;
        h.XLabel.FontSize=12;
        figh=gcf;
        set(0, 'DefaultFigureRenderer', 'painters');
        set(figh,'WindowStyle','normal')
        set(figh, 'PaperUnits', 'inches','Renderer','painters');
        set(figh,'Units','inches','Position',[1 1 2 2])

        exportgraphics(figh,['PinkNoise_fLIF_eta' suffix{a} '.eps'])
        set(figh,'WindowStyle','docked')
    end
    input('r')
end

mdV=diff(slopeV,[],2);
%% Figure 4D (but in the AI files is 5)
clf
%yyaxis left
errorbar(slopeV(:,1),mdV,sertVM,'color',[0 0 0])
ylim([0 1])
ylabel(['\Delta\eta_{LIF}=\beta-\eta_{LIFout}'])
% yyaxis right
% plot(slopeV(:,1),whiteIndex(:,2),'o--r')
% hold on
% plot(slopeV(:,1),whiteIndex(:,1),'o-r')
ylim([0 0.5])
xlim([-1 0])
xlabel('Pink noise exponent ')
%ylabel('White index')
set(gca,'clipping','on')
box off
h=gca;
h.FontSize=10;
h.YLabel.FontSize=12;
h.XLabel.FontSize=12;

figh=gcf;
set(0, 'DefaultFigureRenderer', 'painters');
set(figh,'WindowStyle','normal')
set(figh, 'PaperUnits', 'inches','Renderer','painters');
set(figh,'Units','inches','Position',[1 1 2 2])
exportgraphics(figh,'PinkNoise_fLIF_WhiteIndexSummary.eps')
set(figh,'WindowStyle','docked')





