classdef fracCir
    properties
%         dirdat='/home/fidel/Patricia/Data';
        dirdat='/Users/patriciavazquez/Desktop/Doctorado/6to Semestre/Nuevo/SpAnalysis';
        genInfo='Analysis of fractional order electric circuit'
        SeriesName=[];
        Exp=[];
    end
    methods
        function obj = fracCir(varargin) %the constructor
            %you can send fracCir(dirdat,'/home/fidel/Data')
            if ~isempty(varargin)
                obj.(varargin{1})=varargin{2};
            end
        end
        function obj =loadData(obj,Sname,varargin)
            soption='ONE';
            range=[];
            if ~isempty(varargin)
                soption=upper(varargin{1});
                range=varargin{2};
            end
            
            switch soption
                case 'ONE'
                    dummy=load([obj.dirdat '/' Sname]);
                    obj.SeriesName=Sname;
                    for a=1:length(dummy.DataSt)
                        obj.Exp(a).Data=dummy.DataSt(a).Data;
                    end
                case 'SERIES'
                    c=1;
                    for a=range
                        
                        suffix=[ num2str(a)];
                        
                        dummy=load([obj.dirdat '/' suffix '.mat']);
                        obj.SeriesName=Sname;
                        obj.Exp(c).Data=dummy.DataSt.Data;
                        c=c+1;
                    end
            end
        end
        function obj=AnalizeTrain(obj,varargin)
            %simple read files and analyze automatically
            [obj.dirdat '/' obj.SeriesName]
            dSt=obj.Exp;
            if ~isempty(varargin)
                exp2a=varargin{1};
            else
                exp2a=1:length(dSt);
            end
            
            for a=exp2a
                
                dummy=skpanalisys(dSt(a),0.01,-1.6);
                obj.Exp(a).Analysis=dummy.Analysis;
                
            end
        end
       function obj=AnalizeTrainNW(obj,varargin)
            %simple read files and analyze automatically
            [obj.dirdat '/' obj.SeriesName]
            dSt=obj.Exp;
            if ~isempty(varargin)
                exp2a=varargin{1};
            else
                exp2a=1:length(dSt);
            end
            
            for a=exp2a
                
                dummy=skpanalisys(dSt(a),0.01,3.5);
                obj.Exp(a).Analysis=dummy.Analysis;
                
                dSt(a).Data.Vout=dSt(a).Data.Vout2;
                dummy2=skpanalisys(dSt(a),0.01,3.5);
                obj.Exp(a).Analysis2=dummy2.Analysis;
                
            end
        end
        
        
        function outM1=plotMeanFR(obj,varargin)
            h=gca;
            st2a=1;
            exp2a=[];
            if ~isempty(varargin)
                for vv=1:2:length(varargin)
                    option=upper(varargin{vv});
                    switch option
                        case 'STRUCT'
                            st2a=varargin{vv+1};
                        case 'EXP'
                            exp2a=varargin{vv+1};
                    end
                end
            end
            if isempty(exp2a)
                dSt=obj(st2a).Exp;
            else
                dSt=obj(st2a).Exp(exp2a);
            end
            
            %             dSt=obj.Exp;
            outM=[];
            for a=1:length(dSt)
                thisD=dSt(a).Data;
                thisA=dSt(a).Analysis;
%                 outM=[outM; [a mean(thisD.Vin).*thisD.CV mean(thisA.isifrInst(:,3)) std(thisA.isifrInst(:,3))]];
                outM=[outM; [a mean(thisD.Vin).*0.1 mean(thisA.isifrInst(:,3)) std(thisA.isifrInst(:,3))]];
            
            end
            outM1=sortrows(outM,1);
            
            hold on
            errorbar(outM1(:,2), outM1(:,3),outM1(:,4),'LineWidth',3)
%             text(outM1(:,2),outM1(:,3),[num2str(outM1(:,1))],'FontWeight','bold','FontSize',8);
            title('Mean Firing Rate','FontWeight','bold','FontSize',16)
            xlabel('Current (mA)','FontWeight','bold','FontSize',14)
            ylabel('Firing Rate (Hz)','FontWeight','bold','FontSize',14)
            grid on
            axis('tight')
            ax = gca; % current axes
            ax.FontWeight = 'Bold';
            box off
        end
        function outM=plotAvalanche(obj,varargin)
            h=gca;
            
          
            st2a=1;
            exp2a=[];
            if ~isempty(varargin)
                for vv=1:2:length(varargin)
                    option=upper(varargin{vv});
                    switch option
                        case 'STRUCT'
                            st2a=varargin{vv+1};
                        case 'EXP'
                            exp2a=varargin{vv+1};
                    end
                end
            end
            if isempty(exp2a)
                dSt=obj(st2a).Exp;
            else
                dSt=obj(st2a).Exp(exp2a);
            end
            
            outM=[];
            c=1;
            for a=1:length(dSt)
                thisD=dSt(a).Data;
                thisA=dSt(a).Analysis.Av;
                if ~isempty(thisA.Hist)
                    y2=thisA.Hist(1:18,2);
                    x2=thisA.Hist(1:18,1);
%                     f0=thisA.LogFit;
%                     cval=coeffvalues(f0);
                    
                    %                     plot(log10(x2),log10(y2),'Marker','*','MarkerSize',3,'LineWidth',2,'Color',[1 0 0])
                    %                     hold on
                    %                     plot(log10(x2),f0(log10(x2)),'k','LineWidth',2, 'LineStyle','-.')
                    %                     text(log10(x2(end)),f0(log10(x2(end))),num2str(a));
                    %                     %text(log10(x2(end)),f0(log10(x2(end))),num2str(1000*Alltemp(a).Data.CV*mean(Alltemp(a).Data.Vin)))
                    f=fittype('a*x+b');
                    g=fitoptions(f);
                    p1=8;
                    f0=fit(log10(x2(1:end-p1)),log10(y2(1:end-p1)),f,g);
                    out(a).f0=f0;
                    cval=coeffvalues(f0);
               
                    ax = gca; % current axes
                    ax.FontWeight = 'Bold';
    
                    ax = gca; % current axes
                    ax.FontWeight = 'Bold';
                    
                    %                     figure(12)
                    loglog(x2(1:end-p1),y2(1:end-p1));
                    hold on
                    loglog(x2(1:end-p1),10.^f0(log10(x2(1:end-p1))),'k','LineWidth',2, 'LineStyle','-.')
                    text(x2(end-p1),10.^f0(log10(x2(end-p1))),[num2str(a) '    ' num2str(cval(1))],'FontWeight','bold','FontSize',8);
                    drawnow
                    outM(a).m=[x2,y2];
                    outM(a).slope=cval(1);
                    c=c+1;
                end
            end
            title('Avalanches','FontWeight','bold','FontSize',12)
            xlabel('Spikes','FontWeight','bold','FontSize',10)
            ylabel('Counts','FontWeight','bold','FontSize',10)
            grid
            
%             axis('tight')
            box off
        end
        function outM=plotAvalanche2(obj,varargin)
            h=gca;
            
          
            st2a=1;
            exp2a=[];
            if ~isempty(varargin)
                for vv=1:2:length(varargin)
                    option=upper(varargin{vv});
                    switch option
                        case 'STRUCT'
                            st2a=varargin{vv+1};
                        case 'EXP'
                            exp2a=varargin{vv+1};
                    end
                end
            end
            if isempty(exp2a)
                dSt=obj(st2a).Exp;
            else
                dSt=obj(st2a).Exp(exp2a);
            end
            
            outM=[];
            c=1;
            for a=1:length(dSt)
                thisD=dSt(a).Data;
                thisA=dSt(a).Analysis2.Av;
                if ~isempty(thisA.Hist)
                    y2=thisA.Hist(1:end-2,2);
                    x2=thisA.Hist(1:end-2,1);
                    f0=thisA.LogFit;
                    cval=coeffvalues(f0);
                    
                    %                     plot(log10(x2),log10(y2),'Marker','*','MarkerSize',3,'LineWidth',2,'Color',[1 0 0])
                    %                     hold on
                    %                     plot(log10(x2),f0(log10(x2)),'k','LineWidth',2, 'LineStyle','-.')
                    %                     text(log10(x2(end)),f0(log10(x2(end))),num2str(a));
                    %                     %text(log10(x2(end)),f0(log10(x2(end))),num2str(1000*Alltemp(a).Data.CV*mean(Alltemp(a).Data.Vin)))
                    
                    ax = gca; % current axes
                    ax.FontWeight = 'Bold';
                    
                    %                     figure(12)
                    loglog(x2,y2);
                    hold on
                    loglog(x2,10.^f0(log10(x2)),'k','LineWidth',2, 'LineStyle','-.')
                    text(x2(end),10.^f0(log10(x2(end))),[num2str(a) '    ' num2str(cval(1))],'FontWeight','bold','FontSize',8);
                    drawnow
                    outM(a).slope=cval(1);
                    outM(a).m=[x2,y2];
                    c=c+1;
                end
            end
            title('Avalanches','FontWeight','bold','FontSize',12)
            xlabel('Spikes','FontWeight','bold','FontSize',10)
            ylabel('Counts','FontWeight','bold','FontSize',10)
            grid
            
            axis('tight')
            box off
        end
        
        
        function outM=plotIFRPowerSpectrum(obj,varargin)
            h=gca;
            h.FontWeight = 'Bold';
            h.FontSize = 12;
            st2a=1;
            exp2a=[];
            if ~isempty(varargin)
                for vv=1:2:length(varargin)
                    option=upper(varargin{vv});
                    switch option
                        case 'STRUCT'
                            st2a=varargin{vv+1};
                        case 'EXP'
                            exp2a=varargin{vv+1};
                    end
                end
            end
            if isempty(exp2a)
                dSt=obj(st2a).Exp;
            else
                dSt=obj(st2a).Exp(exp2a);
            end
            
            t=dSt(1).Data.t;
            dt=diff(t(1:2));
            Fs=1/dt;
            %             dSt=obj.Exp(n);
            outM=[];
            f=fittype('a*x+b');
            g=fitoptions(f);
            c=1;
            for a=1:length(dSt)
                thisD=dSt(a).Data;
                isifrInst=dSt(a).Analysis.isifrInst;
                if ~isempty(isifrInst)
                    
                 
%                     dummysp=sparse(dSt(a).Analysis.allSpIndex,1,1,length(t),1);
%                     dummyspC=sparse(dSt(a).Analysis.scSp.cSp.SpIndex,1,1,length(t),1);
%                     dummyspI=sparse(dSt(a).Analysis.scSp.ISp.SpIndex,1,1,length(t),1);
%                     sp0=[];
% 
%                     
%                     
%                     sp0=full(dummysp);
%                     sp0C=full(dummyspC);
%                     sp0I=full(dummyspI);
                    

                    ISI=dSt(a).Analysis.isifrInterp(1:end-2,2);%just the firing rate
                    w1=fft(ISI);
                    w2=w1(1:length(w1)/2);
                    fv=linspace(0,5000,length(t)/2);
                    w2log=log10(abs(w2));
                    hold on
                    plot(log10(fv(1:end-1)),w2log-max(w2log(2:end)),'b','LineWidth',3)
                    
                    Vin=dSt(a).Data.Vin;
                    w1=fft(Vin);
                    w2=w1(1:length(w1)/2);
                    fv=linspace(0,5000,length(t)/2);
                    w2log=log10(abs(w2));
                    hold on
                    plot(log10(fv(1:end)),w2log-max(w2log(2:end)),'r','LineWidth',3)


                    bin=40;
                    p0=1;
                    p2=20;
                    fullifr=dSt(a).Analysis.isifrInterp(1:end,2);%just the firing rate
%                     fullifr=sp0C;
                    fullifr=round(fullifr*1000)./1000;%round to ms.
                    P1=obj.getFFT(fullifr(Fs:end-Fs),Fs,bin);
                    Y=P1.BinData(:,2);
                    fv=P1.BinData(:,1);
                    g.StartPoint=[diff(Y([end 1]))/diff(fv([end 1])) 0.83];
                    fv2=fv(logical(~(Y<=0)));
                    Y2=Y(logical(~(Y<=0)));

                    
                    ax = gca; % current axes
                    ax.FontWeight = 'Bold';
                    Y2=Y2(1:end);
                    fv2=fv2(1:end);
                    plot(log10(fv2),log10(Y2./Y2(2)),'Color','magenta','Marker','*','MarkerSize',3,'LineWidth',3);
                    hold on
                    Y2=Y2(p0:p2);
             
                    fv3=fv2(p0:p2);
                    
                    f0=fit([-2; log10(fv3)],[0;log10(Y2./Y2(2))],f,g);
                    cval=coeffvalues(f0);
                    Y3=10.^f0(log10(fv2));
                    plot(log10(fv2),log10(Y3),'k','LineWidth',2, 'LineStyle','-.')
                           
                    %%%%%%%%%%
                    fv4=log10(fv2);
                    WhiteS=log10(Y3);
                    out.WindexFR=trapz(fv4,1+WhiteS)/3
                    plot(fv4,1+WhiteS)
                    
                    
                    %10.^f0(log10(fv(1)))
                    text(log10(fv(1)),log10(Y3(1))+log10(Y2(2)),[num2str(a) '    ' num2str(cval(1))],'FontWeight','bold','FontSize',10);
                    legend('ISI','Vin','ISI binned')
                    
                    outM(a).spectrumIFR = P1.BinData;
                    outM(a).LogFitIFR=f0;
                else
                    outM(a).spectrumIFR =[];
                    outM(a).LogFitIFR=[];
                end
                Pin=obj.getFFT(thisD.Vin,Fs,10);
                Y=Pin.BinData(:,2);
                fv=Pin.BinData(:,1);
                g.StartPoint=[diff(Y([end 1]))/diff(fv([end 1])) 0.83];
                fv2=fv(logical(~(Y<=0)));
                Y2=Y(logical(~(Y<=0)));
                outM(a).spectrumVin=Pin.BinData;
                outM(a).LogFitVin=[];%f1;
                
                fv=linspace(0,5000,length(t)/2);

                pinknoise=w2log-max(w2log(2:end));
                
                
                %%%%%%%%
                fv=fv(pinknoise>-1);
                pinknoise=pinknoise(pinknoise>-1);
                plot(log10(fv),1+pinknoise,'--g','LineWidth',3)
                out.WindexPN=trapz(log10(fv(2:end)),1+pinknoise(2:end))/3
    
         

                
            end
            %axis('tight')
            hold on
            title('IFR Power Spectrum','FontWeight','bold','FontSize',12)
            xlabel('Frequency (Hz)','FontWeight','bold','FontSize',10)
            ylabel('Dispersion','FontWeight','bold','FontSize',10)
            grid on
            xlim([-2,0])
            ylim([-2,1])
            box off
        end
        
        function outM=plotIFRPowerSpectrum2(obj,varargin)
            h=gca;
            h.FontWeight = 'Bold';
            h.FontSize = 12;
            st2a=1;
            exp2a=[];
            if ~isempty(varargin)
                for vv=1:2:length(varargin)
                    option=upper(varargin{vv});
                    switch option
                        case 'STRUCT'
                            st2a=varargin{vv+1};
                        case 'EXP'
                            exp2a=varargin{vv+1};
                    end
                end
            end
            if isempty(exp2a)
                dSt=obj(st2a).Exp;
            else
                dSt=obj(st2a).Exp(exp2a);
            end
            
            t=dSt(1).Data.t;
            dt=diff(t(1:2));
            Fs=1/dt;
            %             dSt=obj.Exp(n);
            outM=[];
            f=fittype('a*x+b');
            g=fitoptions(f);
            c=1;
            for a=1:length(dSt)
                thisD=dSt(a).Data;
                isifrInst=dSt(a).Analysis2.isifrInst;
                if ~isempty(isifrInst)
                    fullifr=dSt(a).Analysis2.isifrInterp(1:end,2);%just the firing rate
                    fullifr=round(fullifr*1000)./1000;%round to ms.
                    P1=obj.getFFT(fullifr(Fs:end-Fs),Fs,10);
                    Y=P1.BinData(:,2);
                    fv=P1.BinData(:,1);
                    g.StartPoint=[diff(Y([end 1]))/diff(fv([end 1])) 0.83];
                    fv2=fv(logical(~(Y<=0)));
                    fv2b=fv2(10:20);
                    Y2=Y(logical(~(Y<=0)));
                    Y2b=Y2(10:20);
                    fvfitr=(fv2b<1e1).*(fv2b>=2e-1);
                    f0=fit(log10(fv2b(logical(fvfitr))),log10(Y2b(logical(fvfitr))),f,g);
                    cval=coeffvalues(f0);
                    
                    ax = gca; % current axes
                    ax.FontWeight = 'Bold';
                    loglog(fv2,Y2,'Marker','*','MarkerSize',3,'LineWidth',3);
                    hold on
                    loglog(fv,10.^f0(log10(fv)),'k','LineWidth',2, 'LineStyle','-.')
                    text(fv(1),10.^f0(log10(fv(1))),[num2str(a) '    ' num2str(cval(1))],'FontWeight','bold','FontSize',10);
                    
                    
                    outM(a).spectrumIFR = P1.BinData;
                    outM(a).LogFitIFR=f0;
                else
                    outM(a).spectrumIFR =[];
                    outM(a).LogFitIFR=[];
                end
                Pin=obj.getFFT(thisD.Vin,Fs,10);
                Y=Pin.BinData(:,2);
                fv=Pin.BinData(:,1);
                g.StartPoint=[diff(Y([end 1]))/diff(fv([end 1])) 0.83];
                fv2=fv(logical(~(Y<=0)));
                Y2=Y(logical(~(Y<=0)));
                outM(a).spectrumVin=Pin.BinData;
                outM(a).LogFitVin=[];%f1;
            end
            %axis('tight')
            hold on
            title('IFR Power Spectrum','FontWeight','bold','FontSize',12)
            xlabel('Frequency (Hz)','FontWeight','bold','FontSize',10)
            ylabel('Dispersion','FontWeight','bold','FontSize',10)
            grid on
            xlim([10^-1,50])
            ylim([10^-3.5,10^0])
            box off
        end
        
        function out=getSpShapes(obj,a,option,av2plot)
            %option: all, avalanche
            dSt=obj.Exp(a);
            thisD=dSt.Data;
            thisA=dSt.Analysis;
            allV=thisD.Vout;
            t=thisD.t;
            dt=diff(t(1:2));
            %             dt=diff(t(1:2)); %sampling rate
            %             %make a filter
            %             windows=round(0.5e-3/dt);
            %             filterb=(1/windows)*ones(1,windows);
            %             filtera=1;
            %             allVf1=filter(filterb,filtera,allV); %filter data;
            %             allVk1=filter(filterb,filtera,thisD.Vout2);
            out.dt=dt;
            out.V=obj.getSpikesShapes(dSt,'V',option,av2plot);
            out.K=obj.getSpikesShapes(dSt,'K',option,av2plot);
        end
        
        
        function [out] = OscTrace(obj,time,varargin)
            %%% osciloscope depending of bursting
            st2a=1;
            %             time=1/5;
            exp2a=[];
            if ~isempty(varargin)
                for vv=1:2:length(varargin)
                    option=upper(varargin{vv});
                    switch option
                        case 'STRUCT'
                            st2a=varargin{vv+1};
                        case 'EXP'
                            exp2a=varargin{vv+1};
                    end
                end
            end
            if isempty(exp2a)
                exp2a=input(['Introduce Signal: ']);
                dSt=obj(st2a).Exp(exp2a);
            else
                dSt=obj(st2a).Exp(exp2a);
            end
            dt=dSt.Data.dt;
            f=dSt.Data.Frequency;
            freq=floor(1/dt*time);
            isifrInst=dSt.Analysis.isifrInterp(:,1);
            rsISI=reshape(isifrInst(1:freq*round(length(isifrInst)/freq)),[freq round(length(isifrInst)/freq)]);
            Mean=mean(rsISI',2);
            %             disp(['The ISI:' num2str(mean(Mean))])
            %             figure
            %             ax = gca; % current axes
            %             ax.FontWeight = 'Bold';
            t=dSt.Data.t;
            cSpI=dSt.Analysis.scSp.cSp.SpIndex;
            for b=2:length(cSpI)%-2
                a=cSpI(b);
                temp2(b)=mod(t(a)*f,1)*2*pi;
                if temp2(b)>6; temp2(b)=temp2(b)-2*pi;end
                                r2p=(a:a+freq);
                                t2p=t(r2p);
                                cla
                                yyaxis left
                                cla
                                plot(t2p,dSt.Data.Vout(r2p));
                                hold on
                                i2p=cSpI(logical((cSpI>=r2p(1)).*(cSpI<=r2p(end))));
                                plot(t(i2p),dSt.Data.Vout(i2p),'r*')
                                yyaxis right
                                plot(t2p,dSt.Data.Vin(r2p))
                                ylim([min(dSt.Data.Vin(r2p)) max(dSt.Data.Vin(r2p))])
                                drawnow
                                input(['cont @ ' num2str(a) ' ' num2str((temp2(b)))])
            end
            out=a;
            out=mean(temp2(2:end));
        end
        
        function [out] = OscTrace2(obj,varargin)
            %oscioloscope depending of period
            st2a=1;
            
            exp2a=[];
            if ~isempty(varargin)
                for vv=1:2:length(varargin)
                    option=upper(varargin{vv});
                    switch option
                        case 'STRUCT'
                            st2a=varargin{vv+1};
                        case 'EXP'
                            exp2a=varargin{vv+1};
                    end
                end
            end
            if isempty(exp2a)
                exp2a=input(['Introduce Signal: ']);
                dSt=obj(st2a).Exp(exp2a);
            else
                dSt=obj(st2a).Exp(exp2a);
            end
            dt=dSt.Data.dt;
            f=dSt.Data.Frequency;
            time=1/f;
            freq=floor(1/dt*time);
            isifrInst=dSt.Analysis.isifrInterp(:,1);
            rsISI=reshape(isifrInst(1:freq*round(length(isifrInst)/freq)),[freq round(length(isifrInst)/freq)]);
            Mean=mean(rsISI',2);
            disp(['The ISI:' num2str(mean(Mean))])
            figure
            ax = gca; % current axes
            ax.FontWeight = 'Bold';
            %             b=input(['Introduce time: ']);
            t=dSt.Data.t;
            c=1;
            cSpI=dSt.Analysis.scSp.cSp.SpIndex;
            %             for b=1:length(temp)
            %                 a=temp(b);
            trl=round(length(isifrInst)/freq);
            for b=1:trl-1
                a=b;
                %                 if (a-1)*(freq)+1<=cSpI(c) && cSpI(c)<=(a)*freq
                r2p=floor((a-1)*(1/dt*time)+1):floor((a)*1/dt*time);
                t2p=t(r2p);
                cla
                yyaxis left
                cla
                plot(t2p,dSt.Data.Vout(r2p));
                hold on
                i2p=cSpI(logical((cSpI>=r2p(1)).*(cSpI<=r2p(end))));
                plot(t(i2p),dSt.Data.Vout(i2p),'r*')
                yyaxis right
                plot(t2p,dSt.Data.Vin(r2p))
                ylim([min(dSt.Data.Vin(r2p)) max(dSt.Data.Vin(r2p))])
                drawnow
                input(['cont @ ' num2str(a) ' ' num2str((a))])
                %                     c=c+1;
                %                 end
            end
            out=a;
        end
        
        function [out] = OscTrace3(obj,varargin)
            %oscioloscope chirp
            st2a=1;
            
            exp2a=[];
            if ~isempty(varargin)
                for vv=1:2:length(varargin)
                    option=upper(varargin{vv});
                    switch option
                        case 'STRUCT'
                            st2a=varargin{vv+1};
                        case 'EXP'
                            exp2a=varargin{vv+1};
                    end
                end
            end
            if isempty(exp2a)
                exp2a=input(['Introduce Signal: ']);
                dSt=obj(st2a).Exp(exp2a);
            else
                dSt=obj(st2a).Exp(exp2a);
            end
            dt=dSt.Data.dt;
            f=1;
            time=1/f;
            freq=floor(1/dt*time);
            isifrInst=dSt.Analysis.isifrInterp(:,1);
            ddVin=diff(diff(dSt.Data.Vin));
            shiftddVin=circshift(ddVin,-1);
            temp=find(ddVin.*1e10>0); 
            temp2=find(shiftddVin.*1e10<0); 
            index=intersect(temp,temp2);
            figure
            ax = gca; % current axes
            ax.FontWeight = 'Bold';
            t=dSt.Data.t;
            c=1;
            cSpI=dSt.Analysis.scSp.cSp.SpIndex;
            sps=zeros(length(index),30);
            for b=1:length(index)-1
                r2p=index(b):index(b+1);
                t2p=t(index(b):index(b+1)); 
                freq=1/(max(t2p)-min(t2p));
                timeperiod=(t2p-min(t2p))*2*pi*freq;
                i2p=cSpI(logical((cSpI>=r2p(1)).*(cSpI<=r2p(end))));
                timepsp=(t(i2p)-min(t2p))*2*pi*freq;
                sps(b,(1:length(timepsp)))=timepsp;
                cla
                yyaxis left
                cla  
%                 plot(timeperiod,dSt.Data.Vout(r2p));
                hold on
                plot(timeperiod,dSt.Analysis.isifrInterp(r2p,2));
                hold on
                [ft] = sinefit(timeperiod,dSt.Analysis.isifrInterp(r2p,2));
                tSine2=0:0.1:2*pi;
                y=ft(1)+ft(2).*sin(2.*pi.*ft(3).*tSine2+ft(4));
                hold on
                plot(tsine2,y)
%                 plot(timepsp,dSt.Data.Vout(i2p),'r*')
                yyaxis right
                plot(timeperiod,dSt.Data.Vin(r2p))
                ylim([min(dSt.Data.Vin(r2p)) max(dSt.Data.Vin(r2p))])
                drawnow
                timepsp
                input(['cont @ ' num2str(b) ' ' num2str((b))])
                
            end
            out=sps;
        end
        function [out] = DiffusionEntropy(obj,varargin) %S,del,P1,P2)
            st2a=1;
            exp2a=[];
            ISItype=[];
            if ~isempty(varargin)
                for vv=1:2:length(varargin)
                    option=upper(varargin{vv});
                    switch option
                        case 'STRUCT'
                            st2a=varargin{vv+1};
                        case 'EXP'
                            exp2a=varargin{vv+1};
                        case 'ISI'
                            ISItype=varargin{vv+1};
                    end
                end
            end
            if isempty(exp2a)
                dSt=obj(st2a).Exp;
            else
                dSt=obj(st2a).Exp(exp2a);
            end
            for a=1:length(dSt)
                %%%%%%% F Cumulative distribution function
                if ISItype==1
                    isifrInst=dSt(a).Analysis.isifrInst;
                    S=0;
                    S=isifrInst(:,2);
                    
                    P1=1;
                    P2=1;
                elseif ISItype==2
                    isifrInst=dSt(a).Analysis.isifrInterp;
                    S=0;
                    S=isifrInst(1:end,1);
                    P1=1;
                    P2=1;
                    
                elseif ISItype==3
                    isifrInst=dSt(a).Analysis.isifrInst;
                    %                     isifrInst=dSt(a).Data.Vin;
                    S=0;
                    S=isifrInst(1,1:end);
                    P1=1;
                    P2=5;
                else
                    isifrInst=dSt(a).Analysis.isifrInst;
%                     isifrInst=dSt(a).Data.Vin;
                    S=0;
                    S=isifrInst(1:end,2);
%                     
                    
%                     S=isifrInst(:,1);
                    %                     [m,n] = size(S) ;
                    %                     idx = randperm(m) ;
                    %                     b = S ;
                    %                     b(idx,1) = S(:,1);
                    %                     S=b;
                    P1=1;
                    P2=3;
                end

                S1=S-min(S);
                S2=(S1)/max(S1);
                S=S2;
                dispersion=(max(S)-min(S))/20;
                N=length(S); % Using stripes
%                 del=dispersion
                del=0.01;
%                 del=0.15;
                ibeatc=ceil(S./del).*2;
                ibeat=(ibeatc+floor(S./del).*2)/2;
                k=1;
                tau=[];
                ib=[];
                tau(1)=1;
                ib=ibeat;
                for o=2:N
                    if (ib(o)==ib(o-1))
                        tau(k)=tau(k)+1;
                    else
                        k=k+1;
                        tau(k)=1;
                    end
                end
                w=tau;
                %%%Diffusion trajectory- step ahead every time there is an event
                lnt=0;
                den=[];
                den1=[];
                parfor j=1:length(w)
                    den=zeros(w(j)-1,1);
                    den1=[den;1];
                    lnt=[lnt;den1];
                end
                lnt=lnt(2:end);
                X=0;
                X=cumsum(lnt); %Cumulative function distribution
                %%%%%%%%%%%Log Time%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ki=length(X);
                l=floor(log(ki)/log(2))-1; %%%%%
                LogTime=zeros(1,l);
                for i=1:l
                    LogTime(i)=2^i;
                end
                %%%%%%%%%%%%%% Entropy each time
                DE=[];
                DEF=[];
                EEE=[];
                for q=1:length(LogTime)
                    
                    [DEF]= obj.functionEntropy(LogTime(q),X);
                    EEE(:,q)= DEF;
                end
                DE(:,1)=EEE;
                %%%%%%%%%%%% PLOT Difussion Entropy and slope
                x2=0;
                y2=0;
                LogTime=LogTime';
                x2=LogTime(P1:end-P2);
                y2=DE(P1:end-P2);
                y2(logical(y2==0))=1;
                x2(logical(x2==0))=1;
                f=fittype('a*x+b');
                g=fitoptions(f);
                f0=fit(log10(x2),y2,f,g);
                out(a).f0=f0;
                cval=coeffvalues(f0);
                xlog=x2;
                ylog=f0(log10(x2));
                ax = gca; % current axes
                ax.FontWeight = 'Bold';
                semilogx(LogTime,DE,'o')%,'Color',color)
                hold on
                semilogx(x2,f0(log10(x2)),'--k')
                scale=cval(1)/log(10);
                text(xlog(end)+1,ylog(end),[num2str(a), '\delta=', num2str(scale,'%.2f')],'FontWeight','bold','FontSize',8)
                title('Difussion entropy','FontSize',12)
                ylabel('Difussion entropy (S(\tau))','FontWeight','bold','FontSize',10)
                xlabel('Logarithmic time (\tau)','FontWeight','bold','FontSize',10)
                grid on
                out(a).f0=f0;
                out(a).scale=scale;
                out(a).xlog=xlog;
                out(a).ylog=ylog;
            end
        end
        
        
        function [out] = DiffusionEntropy2(obj,varargin) %S,del,P1,P2)
            st2a=1;
            exp2a=[];
            ISItype=[];
            if ~isempty(varargin)
                for vv=1:2:length(varargin)
                    option=upper(varargin{vv});
                    switch option
                        case 'STRUCT'
                            st2a=varargin{vv+1};
                        case 'EXP'
                            exp2a=varargin{vv+1};
                        case 'ISI'
                            ISItype=varargin{vv+1};
                    end
                end
            end
            if isempty(exp2a)
                dSt=obj(st2a).Exp;
            else
                dSt=obj(st2a).Exp(exp2a);
            end
            for a=1:length(dSt)
                %%%%%%% F Cumulative distribution function
                if ISItype==1
                    isifrInst=dSt(a).Analysis2.isifrInst;
                    S=0;
                    S=isifrInst(:,2);
                    P1=1;
                    P2=1;
                elseif ISItype==2
                    isifrInst=dSt(a).Analysis2.isifrInterp;
                    S=0;
                    S=isifrInst(100:end,1);
                    P1=1;
                    P2=1;
                else
                    isifrInst=dSt(a).Analysis2.isifrInst;
                    S=0;
                    S=isifrInst(1:end,2);
                    %                     [m,n] = size(S) ;
                    %                     idx = randperm(m) ;
                    %                     b = S ;
                    %                     b(idx,1) = S(:,1);
                    %                     S=b;
                    P1=1;
                    P2=3;
                end
                dispersion=0.75*std(S);
                N=length(S); % Using stripes
                del=dispersion;
                ibeatc=ceil(S./del).*2;
                ibeat=(ibeatc+floor(S./del).*2)/2;
                k=1;
                tau=[];
                ib=[];
                tau(1)=1;
                ib=ibeat;
                for o=2:N
                    if (ib(o)==ib(o-1))
                        tau(k)=tau(k)+1;
                    else
                        k=k+1;
                        tau(k)=1;
                    end
                end
                w=tau;
                %%%Diffusion trajectory- step ahead every time there is an event
                lnt=0;
                den=[];
                den1=[];
                parfor j=1:length(w)
                    den=zeros(w(j)-1,1);
                    den1=[den;1];
                    lnt=[lnt;den1];
                end
                lnt=lnt(2:end);
                X=0;
                X=cumsum(lnt); %Cumulative function distribution
                %%%%%%%%%%%Log Time%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ki=length(X);
                l=floor(log(ki)/log(2));
                LogTime=zeros(1,l);
                for i=1:l
                    LogTime(i)=2^i;
                end
                %%%%%%%%%%%%%% Entropy each time
                DE=[];
                DEF=[];
                EEE=[];
                for q=1:length(LogTime)
                    
                    [DEF]= obj.functionEntropy(LogTime(q),X);
                    EEE(:,q)= DEF;
                end
                DE(:,1)=EEE;
                %%%%%%%%%%%% PLOT Difussion Entropy and slope
                x2=0;
                y2=0;
                LogTime=LogTime';
                x2=LogTime(P1:end-P2);
                y2=DE(P1:end-P2);
                y2(logical(y2==0))=1;
                x2(logical(x2==0))=1;
                f=fittype('a*x+b');
                g=fitoptions(f);
                f0=fit(log10(x2),y2,f,g);
                out(a).f0=f0;
                cval=coeffvalues(f0);
                xlog=x2;
                ylog=f0(log10(x2));
                ax = gca; % current axes
                ax.FontWeight = 'Bold';
                semilogx(LogTime,DE,'*')%,'Color',color)
                hold on
                semilogx(x2,f0(log10(x2)),'-k')
                scale=cval(1)/log(10);
                text(xlog(end)+1,ylog(end),[num2str(a), '\delta=', num2str(scale,'%.2f')],'FontWeight','bold','FontSize',8)
                title('Difussion entropy','FontSize',12)
                ylabel('Difussion entropy (S(\tau))','FontWeight','bold','FontSize',10)
                xlabel('Logarithmic time (\tau)','FontWeight','bold','FontSize',10)
                grid on
                out(a).f0=f0;
                out(a).scale=scale;
                out(a).xlog=xlog;
                out(a).ylog=ylog;
            end
        end
        
        
        function [out,ifriM,volM]=TestSine(obj,varargin)
            fig=44;
            line='-';
            st2a=1;
            exp2a=[];
            s=2;
            p=1;
            if ~isempty(varargin)
                for vv=1:2:length(varargin)
                    option=upper(varargin{vv});
                    switch option
                        case 'STRUCT'
                            st2a=varargin{vv+1};
                        case 'EXP'
                            exp2a=varargin{vv+1};
                        case 'LINE'
                            line=varargin{vv+1};
                    end
                end
            end
            if isempty(exp2a)
                dSt=obj(st2a).Exp;
            else
                dSt=obj(st2a).Exp(exp2a);
            end
            obj(st2a).SeriesName
            figure(fig)
            subplot(1,s,p)
            ax = gca;
            ax.FontWeight = 'Bold';
            ax.ColorOrderIndex = 1;
            for c=1:length(dSt) %% DATA
                f(c)=dSt(c).Data.Frequency;
                Vin(:,c)=dSt(c).Data.Vin;
                t(:,c)=dSt(c).Data.t;
                Amp=dSt(c).Data.Amplitude;
                dtv=diff(dSt(c).Data.t);
                dt=dtv(1);
                %                 sout3=dSt(c).Analysis;
                %                 IFRi=sout3.isifrInterp(:,2);
                taupt=round((1/f(c)).*(1/dt)); %points in one period
                dvin=floor((1/dt/f(c)+1:1/dt/f(c):size(Vin(:,c)-1/dt/f(c)))'); %periods
                ifriM=[];
                volM=[];
                dummysp=sparse(dSt(c).Analysis.allSpIndex,1,1,length(t),1);
                for a=1:length(dvin)-1
                    ifriM(:,a)=dummysp(dvin(a)+[1:taupt]);
                    volM(:,a)=dSt(c).Data.Vout(dvin(a)+[1:taupt]);
                    %                     ifriM(:,a)=IFRi(dvin(a)+[1:taupt]);
                end
                out(c).ifrM=ifriM(:,1:end);
                out(c).volM=volM(:,1:end);
                sptaub=mean(ifriM(:,1:end),2); %Mean of all IFRi
                sptaubr=sptaub;
                %                 sptaubr=sptaub-mean(sptaub); %Normalice IFRi
                bin1=f(c)/(0.01);
                bin1=round(taupt/400);%.02*10000;%round(taupt/50); %%Binning IFRi.05*5000;.001*10000;%
                mb=floor(length(sptaubr)/bin1);
                mn=bin1*mb;
                b1=squeeze(sum(reshape(sptaubr(1:mn),[bin1 mb]),1));
                sptaubmax=max(b1);
                sptaubmin=min(b1);
                out(c).AmplitudesIFR=(sptaubmax-sptaubmin)/2;
                b2=b1./max(abs(b1));
                b2=b1.*(1/dt)/bin1;
                bint=t(1:bin1:mn,c);
                %%%% Plot
                figure(fig)
                ax.ColorOrderIndex = c;
                ax.FontWeight = 'bold';
                subplot(1,s,p)
                tSine=bint*2*f(c)*pi;
                out(c).b=b2;
                out(c).tSine=tSine;
                h(c)=plot(tSine,b2,line,'LineWidth',2);  %Plot mean of IFRi
                hold on
                [ft] = sinefit(tSine',b2);
                tSine2=0:0.1:2*pi;
                y=ft(1)+ft(2).*sin(2.*pi.*ft(3).*tSine2+ft(4));
                subplot(1,s,p)
                figure(fig)
                
                ax.ColorOrderIndex = c;
                ax.FontWeight = 'bold';
                
                %                 plot(tSine2,y,':')
                xlim([0 2*pi])
                set(gca,'XTick',0:pi/2:2*pi)
                set(gca,'XTickLabel',{'0','\pi/2','\pi','3*\pi/2','2*\pi'},'FontWeight','bold')
                hold on
                xlabel('Period (Radians)','FontWeight','bold','FontSize',10)
                ylabel('IFR Normalized','FontWeight','bold','FontSize',10)
                axis tight
                %%%% Alpha between input and IFRi
                alphaA(c)=ft(4)*2/pi;
                meanIFRAmp=(max(b1)-min(b1))/2;
                AlphaAmp=log(meanIFRAmp/Amp)/log(2*pi*f(c));
                out(c).AlphaAmp=AlphaAmp;
                %%% Sine fractional derivative
                tsin=0:dt:1/f(c);
                alpha=alphaA(c);
                dy=((2*pi*f(c))^(alpha))*Amp*sin(2*pi*f(c)*tsin+(alpha*pi)/2);
                out(c).AmpDer=max(dy); %Derivative amplitude
                out(c).dy=dy;
                out(c).alpha=nonzeros(alphaA(c));
                out(c).f=nonzeros(f(c));
                out(c).maxii=max(b1);
                out(c).minii=min(b1);
                out(c).phase=ft(4);
                out(c).phaseG=ft(4)*360/(2*pi);
                out(c).Amp=(out(c).maxii-out(c).minii)/2;
                out(c).IFR=b2;
                phase(c)=pi/2-ft(4);
            end
            if mean(Vin(:,1))>10; f=f.*1000; end
            figure(fig)
            subplot(1,s,p)
            vin=Vin(1:dvin,c)-mean(Vin(1:dvin,c));
            plot(2*pi*f(c)*t(1:dvin,c),vin./max(vin),'k-','LineWidth',2)
            legend(h,num2str(f'))
            grid on
            %%%% Plot Alpha
            if length(dSt)>1
                figure(fig+1)
                
                subplot(1,s,p)
                ax = gca;
                ax.FontWeight = 'Bold';
                ylabel('Amplitude','FontWeight','bold','FontSize',10)
                ylabel('Phase','FontWeight','bold','FontSize',10)
                plot(f,phase,'*--','LineWidth',3)
                hold on
                grid on
                xlabel('Input frequency (Hz)','FontWeight','bold','FontSize',10)
            end
        end
        
        function [out,ifriM,volM]=TestSine4(obj,varargin)
            fig=44;
            line='-';
            st2a=1;
            exp2a=[];
            s=1;
            p=1;
            if ~isempty(varargin)
                for vv=1:2:length(varargin)
                    option=upper(varargin{vv});
                    switch option
                        case 'STRUCT'
                            st2a=varargin{vv+1};
                        case 'EXP'
                            exp2a=varargin{vv+1};
                        case 'LINE'
                            line=varargin{vv+1};
                    end
                end
            end
            if isempty(exp2a)
                dSt=obj(st2a).Exp;
            else
                dSt=obj(st2a).Exp(exp2a);
            end
            obj(st2a).SeriesName
            figure(fig)
            subplot(1,s,p)
            ax = gca;
            ax.FontWeight = 'Bold';
            ax.ColorOrderIndex = 1;
            for c=1:length(dSt) %% DATA
                f(c)=dSt(c).Data.Frequency;
                Vin(:,c)=dSt(c).Data.Vin;
                t(:,c)=dSt(c).Data.t;
                Amp=dSt(c).Data.Amplitude;
                dtv=diff(dSt(c).Data.t);
                dt=dtv(1);
                
                sout3=dSt(c).Analysis;
                IFRi=sout3.isifrInterp(:,2);
                taupt=round((1/f(c)).*(1/dt)); %points in one period
                dvin=floor((1/dt/f(c)+1:1/dt/f(c):size(Vin(:,c)-1/dt/f(c)))'); %periods
                ifriM=[];
                ifriM2=[];
                volM=[];
                nspikes=[];
%                 dSt2=[];
%                 IFRi=[];
%                 ISIi=[];
                dummysp=sparse(dSt(c).Analysis.allSpIndex,1,1,length(t),1);
                
                for a=1:length(dvin)-1
%                     dSt2.Data.Vin=dSt(c).Data.Vin(dvin(a)+[1:taupt]);
%                     dSt2.Data.Vout=dSt(c).Data.Vout(dvin(a)+[1:taupt]);
%                     dSt2.Data.t=dSt(c).Data.t([1:taupt]);
%                     dummy(a,:)=skpanalisys(dSt2,0.01,4);
%                     IFRi(a,:)=dummy(a).Analysis.isifrInterp(:,2);
%                     ISIi(a,:)=dummy(a).Analysis.isifrInterp(:,1);
                    ifriM(:,a)=dummysp(dvin(a)+[1:taupt]);
                    volM(:,a)=dSt(c).Data.Vout(dvin(a)+[1:taupt]);
                    ifriM2(:,a)=IFRi(dvin(a)+[1:taupt]);
                    nspikes(a)=sum(dummysp(dvin(a)+[1:taupt]),1);
                    std2(a)=std(IFRi(dvin(a)+[1:taupt]));
                    Mean(a)=mean(IFRi(dvin(a)+[1:taupt]));
                    std3(a)= ((sum((IFRi(dvin(a)+[1:taupt])-mean(IFRi)).^2))/(length(IFRi(dvin(a)+[1:taupt]))-1)).^(1/2);
                end
                
%%%%%%%%%%% Option 1
                temp=mean(ifriM,2);
                sptaubr=temp;%(100:end-100);
%                 sptaubr=sptaub-mean(sptaub); %Normalice IFRi
%                 bin1=.01*10000;
                bin1=round(taupt/50);%.02*10000;%round(taupt/50); %%Binning IFRi.05*5000;.001*10000;%
                mb=floor(length(sptaubr)/bin1);
                mn=bin1*mb;
                b1=squeeze(sum(reshape(sptaubr(1:mn),[bin1 mb]),1));
                b2=b1/bin1*10000;
                bint=t(1:bin1:mn,c);
                tSine=bint*2*f(c)*pi;

                figure(10)
                ax=gca;
                ax.ColorOrderIndex = c;
                [ft] = sinefit(tSine',b2);
                y2=ft(1)+ft(2).*sin(2.*pi.*ft(3).*tSine+ft(4));
                hs(c)=plot(tSine,b2,'LineWidth',3)
                hold on
                plot(tSine,y2,'k:','LineWidth',1)
                xlim([0 2*pi])
                set(gca,'XTick',0:pi/2:2*pi)
                set(gca,'XTickLabel',{'0','\pi/2','\pi','3*\pi/2','2*\pi'},'FontWeight','bold')
                legend(hs,num2str(f'))
                xlabel('Period (radians)','FontWeight','bold','FontSize',12)
                ylabel('IFR interpolate','FontWeight','bold','FontSize',12)
                grid on

                
%%%%%%%%%%%%%%Option 2
                temp2=mean(ifriM2,2);
                t2=dSt(c).Data.t([1:taupt]);
                [ft] = sinefit(t2*2*f(c)*pi,temp2);
                tSine2=0:0.1:2*pi;
                y=ft(1)+ft(2).*sin(2.*pi.*ft(3).*tSine2+ft(4));
                A(c)=2*ft(2);
                Period(c)=1/f(c);
                PhaseFreq(c)=ft(3);
                Phase(c)=ft(4);
                if Phase(c)>6; Phase(c)=Phase(c)-2*pi; end
                %%%% Plot
                figure(fig)
                ax.ColorOrderIndex = c;
                subplot(1,s,p)
                out(c).b=b2;
                out(c).tSine=tSine;
                t2=dSt(c).Data.t([1:taupt]);
                h(c)=plot(t2*2*f(c)*pi,temp2,line,'LineWidth',3);
                hold on
                plot(tSine2,y,'k:','LineWidth',1)
            end
            
%%%%%%%%%%%%%Option 1            
            figure(10)
            yyaxis right
            vin=Vin(1:dvin,c);
            plot(2*pi*f(c)*t(1:dvin,c),vin.*.1,'k-','LineWidth',3)
            ylabel('I_{input} (mA)','FontWeight','bold','FontSize',12)
            axis tight
            legend(hs,num2str(1./f'))
            grid on
            
%%%%%%%%%%%%%Option 2            
            figure(fig)
            yyaxis right
            vin=Vin(1:dvin,c);
            plot(2*pi*f(c)*t(1:dvin,c),vin.*.1,'k-','LineWidth',3)
            ylabel('I_{input} (mA)','FontWeight','bold','FontSize',12)
            axis tight
            legend(h,num2str([1./f' f']))
            grid on          
            
%%%%%%%%%%%%%%%% phase and gain
            Phase2=Phase*360/(2*pi);
            pp=[f',Phase2'];
            pp2=sortrows(pp);
            figure(2)
            subplot 121
            hold on
            plot(pp2(:,1),pp2(:,2),'-*','LineWidth',3,'MarkerSize',10)
            ax = gca;
            ax.FontWeight = 'Bold';
            xlabel('Period length (s)','FontWeight','bold','FontSize',12)
            ylabel('Phase lead (deg)','FontWeight','bold','FontSize',12)
            grid on
            AP=[f',A'];
            AP2=sortrows(AP);
            figure(2)
            subplot 122
            hold on
            plot(AP2(:,1),AP2(:,2),'-*','LineWidth',3,'MarkerSize',10)
            ax = gca;
            ax.FontWeight = 'Bold';
            xlabel('Period length (s)','FontWeight','bold','FontSize',10)
            ylabel('Gain (Hz/mA) (peak to peak)','FontWeight','bold','FontSize',10)
            grid on
            
            figure(3)
            hold on
            subplot 121
            loglog(pp2(:,1),pp2(:,2),'-*','LineWidth',3,'MarkerSize',10)
            ax = gca;
            ax.FontWeight = 'Bold';
            xlabel('Period length (s)','FontWeight','bold','FontSize',12)
            ylabel('Phase lead (deg)','FontWeight','bold','FontSize',12)
            
            fv=pp2(:,1);
            Y=pp2(:,2);
            fv2=fv(logical(~(Y<=0)));
            Y2=Y(logical(~(Y<=0)));
            fvfitr=(fv2<10).*(fv2>=0);
            f2=fittype('a*x+b');
            g=fitoptions(f2);
            f0=fit(log10(fv2(logical(fvfitr))),log10(Y2(logical(fvfitr))),f2,g);
            cval=coeffvalues(f0);
            hold on
            loglog(fv2,10.^f0(log10(fv2)),'k','LineWidth',2, 'LineStyle','-.')
            text(fv2(end),10.^f0(log10(fv2(end))),[num2str(cval(1))],'FontWeight','bold','FontSize',14);
            grid on
            figure(3)
            subplot 122
            
            hold on
            loglog(AP2(:,1),AP2(:,2),'-*','LineWidth',3,'MarkerSize',10)
            ax = gca;
            ax.FontWeight = 'Bold';
            xlabel('Period length (s)','FontWeight','bold','FontSize',10)
            ylabel('Gain (Hz/mA) (peak to peak)','FontWeight','bold','FontSize',10)
            grid on
   
        end
    end
    methods(Static)
        function out=getFFT(fullifr,Fs,bs)
            L=length(fullifr);
            relifr=(fullifr-mean(fullifr));
            Y = fft(relifr);
            P2 = abs(Y/L);
            P1 = P2(1:round(L/2)+1);
            P1(2:end-1) = 2*P1(2:end-1);
            fv=(Fs*(0:round(L/2))/L);
            p.lowV=1e-2;%2/(dt*L);%lowest possible frequency
            p.highV=1e0;%max is max(fv)
            p.data(:,1)=fv';
            p.data(:,2)=P1';
            p1bin=LogBinData(p,bs);
            out.BinData=p1bin.BinData;
            out.OriData=p.data;
        end
        
        
        
        
        
        function H=plotPhasePlane(w,x,sorcx,morallx,dorlx,y,sorcy,morally,dorly)
            xtxt=[];
            ytxt=[];
            dt=w.dt;
            xtr=w.(x).(sorcx).(morallx);
            ytr=w.(y).(sorcy).(morally);
            
            
            if strcmp(upper(dorlx),'DIFF') && strcmp(upper(dorly),'DIFF')
                thisTrx=diff(xtr)/dt;
                thisTry=diff(ytr)/dt;
                xtxt=['d' x '/dt'];
                ytxt=['d' y '/dt'];
            elseif strcmp(upper(dorlx),'LIN') && strcmp(upper(dorly),'DIFF')
                thisTrx=xtr(1:end-1,:);
                thisTry=diff(ytr)/dt;
                xtxt=x;
                ytxt=['d' y '/dt'];
            end
            
            
            H=plot(thisTrx,thisTry);
            xlabel(xtxt);ylabel(ytxt);
            box off
        end
        function out=detectSpikes(Stin,thup,ctref)
            Vout=Stin.Data.Vout;
            Vin=Stin.Data.Vin;
            t=Stin.Data.t;
            difft=diff(t);
            dt=difft(1);
            [sppFunction]=sq2delta(Vout',thup); %  detect the spike in the max point
            spp=sppFunction';
            sppf=find(spp); %spike positions
            %timespp=t(sppf); %spike times
            
            if isempty(sppf)
                out=[];
                return
            end
            
            thisSp=sppf(1);
            spindex=[];
            spindex(1)=thisSp;
            %clean for multiple hits on the same spike
            c=2;
            for a=2:length(sppf)
                if (sppf(a)-thisSp)>(ctref/dt) %larger than a refractory period
                    spindex(c)=sppf(a);
                    thisSp=sppf(a);
                    c=c+1;
                end
            end
            spp2=zeros(size(Vout));
            spp2(spindex)=1; %the spike times after cleaning with refractory period
            out=spindex;
        end
        
        function out=getSpikesShapes(dSt,vork,option,av2plot)
            thisD=dSt.Data;
            thisA=dSt.Analysis;
            t=thisD.t;
            dt=diff(t(1:2));
            if upper(vork)=='V'
                allV=thisD.Vout;
            elseif upper(vork)=='K'
                allV=thisD.Vout2;
            end
            windows=round(0.5e-3/dt);
            filterb=(1/windows)*ones(1,windows);
            filtera=1;
            Vf=filter(filterb,filtera,allV); %filter data;
            
            switch upper(option)
                case 'ALL'
                    cSpIndex=dSt.Analysis.scSp.cSp.SpIndex;
                    sSpIndex=dSt.Analysis.scSp.sSp.SpIndex;
                    cSpRange=dSt.Analysis.scSp.cSp.SpRange;
                    sSpRange=dSt.Analysis.scSp.sSp.SpRange;
                case 'AVALANCHE'
                    cSpIndex=dSt.Analysis.Av.cSpIndex(av2plot,:);
                    sSpIndex=dSt.Analysis.Av.sSpIndex{av2plot};
                    allindex=dSt.Analysis.scSp.cSp.SpIndex;
                    drange=dSt.Analysis.scSp.cSp.SpRange;
                    cSpRange=[];
                    for a=cSpIndex
                        cSpRange=[cSpRange; drange(logical(allindex==a),:)];
                    end
                    allindex=dSt.Analysis.scSp.sSp.SpIndex;
                    drange=dSt.Analysis.scSp.sSp.SpRange;
                    sSpRange=[];
                    for a=sSpIndex
                        sSpRange=[sSpRange; drange(logical(allindex==a),:)];
                    end
            end
            %cSpRange=dSt.Analysis.scSp.cSp.SpRange;
            maxcSpL=max(diff(cSpRange,[],2));
            minxcSpL=min(diff(cSpRange,[],2));
            cSpShape=zeros(maxcSpL,length(cSpIndex));
            
            
            %sSpRange=dSt.Analysis.scSp.sSp.SpRange;
            maxsSpL=max(diff(sSpRange,[],2));
            minsSpL=min(diff(sSpRange,[],2));
            sSpShape=zeros(maxsSpL,length(sSpIndex));
            
            for a=1:length(cSpIndex)
                thisR=cSpRange(a,1):cSpRange(a,2);
                thisR([1 end])
                cSpShape(1:length(thisR),a)=Vf(thisR);
            end
            
            for a=1:length(sSpIndex)
                thisR=sSpRange(a,1):sSpRange(a,2);
                sSpShape(1:length(thisR),a)=Vf(thisR);
            end
            out.cSp.Shapes=cSpShape(1:minxcSpL,:);
            out.cSp.Mean=mean(cSpShape(1:minxcSpL,:),2);
            out.sSp.Shapes=sSpShape(1:minsSpL,:);
            out.sSp.Mean=mean(sSpShape(1:minsSpL,:),2);
        end
        
        function out=thicknessSpike(thisV0,th)
            %th=0.5 for spikes
            %calculate thickness of spikes
            bsp1=find((diff(thisV0>th)==1));
            bsp2=find((diff(thisV0>th)==-1));
            out=bsp2(2:end)-bsp1(1:end-1);
        end
        
        function out= functionEntropy(del,tau)
            p=length(tau);
            j=1;
            for i=1:1:(p-del)
                XF(j)=tau(i+del)-tau(i);
                j=j+1;
            end
            XF=XF(XF~=0);
            %%%%%%%%%%Entropy %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nbins = min(XF):1:max(XF);%te/50;
            nbins(nbins==0)=[];
            [counts]=hist(XF,nbins); % ,centers
            
            counts = counts(counts ~= 0); % remove zeros
            
            out=-sum((counts./sum(counts)).*log((counts./sum(counts))));
        end
    end
end
