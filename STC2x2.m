% Uncoded Alamouti's Space Time Code (STC) 
% with 2 Tx and 2 Rx antennas
% Symbol Error Rate (SER) & Bit Error Rate (BER) 
% Monte-Carlo simulation
close all 
clear all
clc
drawnow
% Setup params
SNRs=0:2:20; % SNRs under consideration [dB]
N=2^12; % Block length
channelModels=[ 2]; % 1 - AWGN, 2 - Rayleigh 
ModulationTypes=[2]; % 1 - BPSK, 2 - QPSK, 3- 8PSK, 4 - 16-QAM
%
legends={'STC 2x2'};
tol=1e-12;
MaxRealizations=1e2; 
MaxSymbolErrors=N*MaxRealizations; % Maximal number of symbol errors per scenario
ModulationsNames={'BPSK','QPSK','8PSK','16 QAM'};
ChannelmodelNames={'AWGN','Rayleigh'};
rhos=10.^(-SNRs/20);
%
[CM,MT,RHO]=ndgrid(channelModels,ModulationTypes,rhos);
ScenarionsNum=numel(RHO);
SymErrs=zeros(numel(ModulationTypes),numel(rhos)); % ChannelModelIndex,ModulationTypeIndex
BitErrs=SymErrs;
RealizedSymbolsNum=zeros(size(SymErrs)); % 
for ScenarioIndex=1:ScenarionsNum
    % Progress report
    disp(['Processing scenaio ' num2str(ScenarioIndex) '/' num2str(ScenarionsNum) ' started at ' datestr(now)])
    drawnow
    % Scenario parameters extraction
    rho=RHO(ScenarioIndex);
    [ChannelModelIndex,ModulationTypeIndex,RhoIndex]=ind2sub(size(RHO),ScenarioIndex);
    channelModel=channelModels(ChannelModelIndex);
    ModulationType=ModulationTypes(ModulationTypeIndex);
    switch ModulationType
        case 1 % BPSK dmin=2
            BPS=1; % Bits Per Symbol
            BitMapping=[1 0]';
            Constellation=(-1).^BitMapping;             % Constellation design               
        case 2 % QPSK dmin=1/sqrt(2)
            BPS=2; % Bits Per Symbol
            BitMapping=[0 1 3 2]'; % Grey code            
            Constellation=(-1i).^BitMapping;    
        case 3
            BPS=3; % Bits Per Symbol
            BitMapping=[0 1 3 2 6 7 5 4];
            Constellation=exp(1i*2*pi*(0:7).'/8);
        case 4 % 16-QAM
            BPS=4;
            BitMapping1D=[0 1 3 2]'; % Grey code
            Constellation1D=(-3:2:3)';
            BitMapping=BitMapping1D*ones(1,4)+4*ones(4,1)*BitMapping1D.';
            Constellation=Constellation1D*ones(1,4)+1i*ones(4,1)*Constellation1D.';
            %
            Constellation=Constellation(:);
            BitMapping=BitMapping(:).';
    end
    SymbolsNum=N/2^(BPS-1);
    assert(mod(SymbolsNum,2)==0)
    Es=mean(abs(Constellation).^2);
    Ms=mean(Constellation);    
    if abs(Ms)>tol || abs(sqrt(Es)-1)>tol
        Constellation=(Constellation-Ms)/sqrt(Es); % Normalize constellation   
        Es=mean(abs(Constellation).^2);
        Ms=mean(Constellation);    
        disp(['Constellation normalization check for scenario ' num2str(ScenarioIndex) '  - Var=' num2str(Es) ' mean=' num2str(mean(Constellation))])
    end
    % Pre allocation
    y=zeros(2,SymbolsNum/2);
    MSEsol=zeros(2,SymbolsNum/2);
    % Monte Carlo simulation
    RealizationIndex=0;
    KeepGoing=true;
    while KeepGoing
        RealizationIndex=RealizationIndex+1;
        % ******
        % * Tx *
        % ******        
        TxSymbolIndices=randi(2^BPS,[SymbolsNum,1]);
        TxBits=BitMapping(TxSymbolIndices);
        TxBits=int2bit(TxBits',BPS); 
        TxBits=TxBits(:);
        s=Constellation(TxSymbolIndices);
        % *******
        % * PHY *
        % *******
        s2=reshape(s,2,SymbolsNum/2); % Reorder to fit STC
        n=(randn(4,SymbolsNum/2)+1i*randn(4,SymbolsNum/2))/sqrt(2); % noise
        switch channelModel
            case 1 % AWGN
                h0=ones(2,SymbolsNum/2);
                h1=ones(2,SymbolsNum/2);
            case 2 % Rayleigh
                h0=(randn(2,SymbolsNum/2)+1i*randn(2,SymbolsNum/2))/sqrt(2) ; % Flat fading                
                h1=(randn(2,SymbolsNum/2)+1i*randn(2,SymbolsNum/2))/sqrt(2) ; % Flat fading                
            otherwise
                error("invalid channel model.")
        end 
        y0(1,:)=([1 1]/sqrt(2))*(s2.*h0)+rho*n(1,:);
        y0(2,:)=([1 -1]/sqrt(2))*((s2).*conj(flipud(h0)))+rho* (n(2,:));
        y1(1,:)=([1 1]/sqrt(2))*(s2.*h1)+rho*n(3,:);
        y1(2,:)=([1 -1]/sqrt(2))*((s2).*conj(flipud(h1)))+rho* (n(4,:));
        % ******
        % * Rx * 
        % ******
        h0abs2=sum(abs(h0).^2,1);%/sqrt(2);
        h1abs2=sum(abs(h1).^2,1);%/sqrt(2);
        habs2=h0abs2+h1abs2;
        MSEsol(1,:)=sum([conj(h0(1,:)); h0(2,:);conj(h1(1,:)); h1(2,:) ].*[y0;y1],1)./habs2;
        MSEsol(2,:)=sum([conj(h0(2,:)); -h0(1,:) ; conj(h1(2,:)); -h1(1,:)].*[y0;y1],1)./habs2;          
        MSEsol2=MSEsol(:).';
        [MSEval,MLindex]=min(abs( ones(size(Constellation))*(MSEsol2) - Constellation*ones(1,SymbolsNum) ).^2,[],1);
        MLindex=MLindex';
        RxBits=BitMapping(MLindex);
        RxBits=int2bit(RxBits',BPS); 
        RxBits=RxBits(:);
        % ********
        % * EVAL *
        % ********
        % Evaluate preformance
        SymErrs(ModulationTypeIndex,RhoIndex)=sum(MLindex~=TxSymbolIndices);
        BitErrs(ModulationTypeIndex,RhoIndex)=sum(RxBits~=TxBits);
        % ******
        % * MC *
        % ******
        KeepGoing=RealizationIndex<MaxRealizations && (SymErrs(ModulationTypeIndex,RhoIndex)<MaxSymbolErrors);
    end
    RealizedSymbolsNum(ModulationTypeIndex,RhoIndex)=RealizationIndex*SymbolsNum;
end
% Analysis
SER=SymErrs./RealizedSymbolsNum;
BER=BitErrs./RealizedSymbolsNum;
% 
for ModulationTypeIndex=1:numel(ModulationTypes)
    ModulationType=ModulationTypes(ModulationTypeIndex );
    hSym=figure('Name',['Symbols ' ModulationsNames{ModulationType} ]);
    hBit=figure('Name',['Bits ' ModulationsNames{ModulationType}] );
    channelModel=channelModels(ChannelModelIndex);
    LineStyle=num2linestyle(1);
    Marker=num2marker(1);
    %
    figure(hSym);        
    semilogy(SNRs, (squeeze(SER(ModulationTypeIndex,:))),strcat(LineStyle, Marker) )
    hold on
    figure(hBit)
    semilogy(SNRs, (squeeze(BER(ModulationTypeIndex,:))),strcat(LineStyle, Marker) )
    hold on
    figure(hSym);  
    legend(legends{:})
    title(['STC: SER vs SNR for '  ModulationsNames{ModulationType} ' modulation.'  ])
    hold off
    grid on
    xlabel('SNR [dB]')
    ylabel('SER')
    figure(hBit);  
    legend(legends{:})
    title(['STC: BER vs SNR for '  ModulationsNames{ModulationType} ' modulation.'  ])
    hold off
    grid on
    xlabel('SNR [dB]')
    ylabel('BER')
end
