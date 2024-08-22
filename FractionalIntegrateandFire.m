
%% this funcion is used in fractional neuron integration. It integrates the fractional derivative and  the  voltage v at each time t.

function out=FractionalIntegrateandFire(NetProp,SimProp)
%rng(SimProp.rseed)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rng(SimProp.rseed)
winsize=SimProp.WinSize;
Ncells=NetProp.Ncells;
alph=SimProp.alpha;%-0.5*rand([Ncells 1]);
refrac=NetProp.iRefrac/SimProp.dt;%steps
vth=NetProp.vTh;
vpeak=NetProp.vPeak;
v0=SimProp.v0;
vrest=SimProp.vrest;
% Delays=SimProp.iDelays;
%Here are the tuning parameters
sig=SimProp.Namp;
%gam=SimProp.gam;
%kc=SimProp.kc;

rm=NetProp.Rm;
taum=NetProp.TauM;

dt=SimProp.dt;
weights=SimProp.iWeights;%.*sign(rand(Ncells,Ncells)-0.5).*~eye(Ncells,Ncells);
% tausyn=NetProp.TauSyn;
 aux2=sign(2*(rand(Ncells,length(SimProp.t),1)-0.5));
    %This code is the part of (S+sig*ep)
    Iinj=SimProp.IinjDC';
    %Iinj=(SimProp.IinjDC'.*ones(Ncells,length(SimProp.t),1)'+sig*aux2');
%     Iinj(1:3*winsize)=0;
    %Iinj(round(length(Iinj)/2)+[0:5000])=0;
    %Iinj(round(length(Iinj)/2):end)=Iinj(1)/3;
t=SimProp.t;
v=zeros(length(t),Ncells);
v(1,:) = v0(1,:);
PSP=v(1,:);
sp=zeros(length(t),Ncells);
sp(1,:)=0*(rand([1,Ncells])<0.3); %initialize spikes

spnew=sp;
PSPnew=PSP;
vnew=v(1,:);
vold=v(1,:);
out.Voltage(1,:)=vnew(1,:);
out.Spikes(1,:)=spnew(1,:);
out.Time(1)=t(1);
out.PSP(1,:)=PSPnew(1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%isstillrefrac=zeros(1,Ncells);%1 if inside refrac period
rfcounter=zeros(1,Ncells);
connectM=~(weights==0);
% setup fractional weights


%Memory weights Caputo based on winsize

NN=round(winsize./dt);%length(winsize);
NN1=NN-1;
nn=1:NN-1;
WCoet=(NN+1-nn).^(1-alph)-(NN-nn).^(1-alph);

%The memory weights Caputo for the entire simulation.
NNfull=length(t);
nnfull=1:NNfull-1;
%WCoet=(NNfull+1-nnfull).^(1-alph)-(NNfull-nnfull).^(1-alph);
% Inhib = zeros(length(t),Ncells);

%Memory weights GL method
c1=ones(NN,1)'; 
%c1=cumprod([alph, (1-(alph+1)./(2:(winsize)))]);
  c2(1)=alph; 
  for n=2:NN
      c2(n)=(1-(1+alph)/n)*c2(n-1);
  end
%flip the weigths to match the integration vector
GLweights=(c2);

%kr = dt^alph*gamma(2-alph);     %  the kernel   from the fractional derivative and  weighted  the markovian term
dtalpha=dt^alph;
dtalphagamma=(dt^alph)*gamma(2-alph);
LengthT=length(t);

myweights=connectM.*weights;

intAlg=SimProp.Algorithm;%'GL';

switch upper(intAlg)
    case 'GL'
        WCoeF=@(wstart,wend)GLweights(wend-wstart:wend);
        DeltaHF=@(d1start,d1end,d2start,d2end,vin)vin(d1end-d1start:d1end,:);
        intF=@(vold,Ihere,wcoe,deltah)dtalpha*(-(vold-vrest)+rm*Ihere)/taum+wcoe*deltah;
    case 'CAPUTO'
        WCoeF=@(wstart,wend)WCoet(wend-wstart:wend);
        DeltaHF=@(d1start,d1end,d2start,d2end,v)v(d1end-d1start:d1end,:)-v(d2end-d2start:d2end,:);
        corrfacF=@(wfrac,a,v,wsize)v(a-winsize,:).*wfrac;
        intF=@(vold,Ihere,wcoe,deltah)dtalphagamma*(-(vold-vrest)+rm*Ihere)/taum + vold - wcoe*deltah;
        intFcorr=@(vold,Ihere,wcoe,deltah,corrfac)dtalphagamma*(-vold+rm*Ihere)/taum + vold - (wcoe*deltah+corrfac);
end
WCoewinsize=WCoeF(winsize-1,winsize);
%the leftout weights when using a window
%WCoeoutside=WCoeF(LengthT-winsize-1,LengthT-winsize);
%WCoeall=WCoeF(LengthT-1,LengthT);
%wcoutfrac=sum(WCoeoutside)./sum(WCoeall);
preT=1:2;
vold=v(1,:);
% for a=preT
%     aparriving=sp(a,:)*(myweights);
%     PSPnew=aparriving;
%     Ihere=(PSPnew+Iinj(a,:));
%     vdummy=dt*(-(vold-vrest)+rm*Ihere)/taum+vold;
% 
%     vold=v(a+1,:);
% end
fii=zeros(size(v));
fii(1:2)=Iinj(1:2);
Delta=diff(v(1:preT(end)+1,:));
for a=1:length(t)-1
    aparriving=sp(a,:)*(myweights);
    PSPnew=aparriving;
    Ihere=(PSPnew+Iinj(a,:));
    if alph<1 %(a>preT(end)) && 
        if a>(NN+1)
            %Delta=[Delta(2:NN1); diff(v(a-2:a-1))];
            %Memory=WCoet*Delta;
            MemoryGL=c2*v(a:-1:a-NN+1);

        elseif a<=NN
            %Memory=WCoet(end-a+preT(end)+1:end)*diff(v(1:a-1));
            MemoryGL=c2(1:a)*(v(a:-1:1));

            %MemoryGL=c2(1:a)*(fii(a:-1:1));
        elseif a==(NN+1)
            %Delta =diff(v(a-NN:a-1));
            %Memory=WCoet*Delta;
            MemoryGL=c2*v(a:-1:2);

            %MemoryGL=c2(1:a)*(fii(a:-1:1));
        end
        %vdummy=dtalphagamma*(-(vold-vrest)+(rm).*Ihere)./(taum) + vold - Memory
        %fii(a+1)=(dtalpha*Ihere + MemoryGL);
        %vdummy=dt*(-(vold-vrest)+rm*fii(a).*(a*dt).^-0)/taum+vold;
        vdummy=dtalpha*(-(vold-vrest)+(rm).*Ihere)./(taum) + MemoryGL;
    else
        vdummy=dt*(-(vold-vrest)+rm*Ihere)/taum+vold;
    end

    %test=(((vdummy>vth) + (isstillrefrac==1))>0) ;%v larger than vth or inside refrac period
    test=(vdummy>vth).*(rfcounter==0);
    sp(a+1,:)=test;
    v(a+1,logical(test))=vpeak;%add a spike
    v(a+1,logical((~test).*(~rfcounter)))=vdummy(logical((~test).*(~rfcounter)));
    v(a+1,logical((rfcounter>0)))=vrest;%vth.*(1-rfcounter(logical(rfcounter>0))./refrac);%this is a hack!
    rfcounter=rfcounter+test+(rfcounter>0);
    rfcounter=rfcounter.*(rfcounter<refrac);%if you pass the counter limit, reset it.
   if ~rem(a+1,round(20000/dt))% if round(1/dt) >trefrac you might loose spikes
       display(['another sec ' num2str(a*dt)])
   end
    vold=v(a+1,:);
end   

out.Voltage=(v);
out.Spikes=(sp);
% out.Iinj=Iinj(1:10:end-1,:);
out.Iinj=Iinj;

