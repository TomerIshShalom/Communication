% Uncoded Eigen beamforming (with one Rx antenna)
% Symbol Error Rate (SER) & Bit Error Rate (BER) 
% Monte-Carlo simulation
close all
clear all
clc
drawnow
% Setup params
SNRs=0:20; % SNRs under consideration [dB]
N=2^12; % Block length
channelModels=[ 2]; % 1 - AWGN, 2 - Rayleigh 
ModulationTypes=[2]; % 1 - BPSK, 2 - QPSK, 3- 8PSK, 4 - 16-QAM
TxNs=[ 4]; % #Rx antennas
%
iif  = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();
RXNnames = arrayfun(@(x) iif(x==1, 'SISO',true ,['MRC 1x' num2str(x) ]),TxNs,'UniformOutput',false);
tol=1e-12;
MaxRealizations=1e5; 
MaxSymbolErrors=MaxRealizations*N; % Max symbols error events number
ModulationsNames={'BPSK','QPSK','8PSK','16 QAM'};
ChannelmodelNames={'AWGN','Rayleigh'};
rhos=10.^(-SNRs/20);
%
[CM,MT,RHO,TXN]=ndgrid(channelModels,ModulationTypes,rhos,TxNs);
ScenarionsNum=numel(RHO);
SymErrs=zeros(numel(TxNs),numel(ModulationTypes),numel(rhos)); % ChannelModelIndex,ModulationTypeIndex
BitErrs=SymErrs;
RealizedSymbolsNum=zeros(size(SymErrs)); % 
for ScenarioIndex=1:ScenarionsNum
    % Progress report
    disp(['Processing scenaio ' num2str(ScenarioIndex) '/' num2str(ScenarionsNum) ' started at ' datestr(now)])
    drawnow
    % Scenario parameters extraction
    rho=RHO(ScenarioIndex);
    [ChannelModelIndex,ModulationTypeIndex,RhoIndex,TxNindex]=ind2sub(size(RHO),ScenarioIndex);
    channelModel=channelModels(ChannelModelIndex);
    ModulationType=ModulationTypes(ModulationTypeIndex);
    TxN=TxNs(TxNindex); % #Rx antennas
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
    Es=mean(abs(Constellation).^2);
    Ms=mean(Constellation);    
    % Verify constellation statistics
    if abs(Ms)>tol || abs(sqrt(Es)-1)>tol
        Constellation=(Constellation-Ms)/sqrt(Es); % Normalize constellation
        Es=sum(abs(Constellation).^2);
        Ms=mean(Constellation);
        disp(['Constellation normalization check for scenario ' num2str(ScenarioIndex) '  - Var=' num2str(Es) ' mean=' num2str(mean(Constellation))])
    end
    %
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
        s=ones(TxN,1)*Constellation(TxSymbolIndices).';
        % Feedback 
        switch channelModel
            case 1 % AWGN
                h=ones(TxN,SymbolsNum);
            case 2 % Rayleigh
                h=(randn(TxN,SymbolsNum)+1i*randn(TxN,SymbolsNum))/sqrt(2) ; % Flat fading
            otherwise
                error("invalid channel model.")
        end   
        % Precoding  
        habs=sqrt( sum(abs(h).^2,1) );
        w=conj(h)./( ones(TxN,1)* habs );
        ps=w.*s;
        % *******
        % * PHY *
        % *******
        n=(randn(1,SymbolsNum)+1i*randn(1,SymbolsNum))/sqrt(2); % noise    
        y=sum(h.*ps,1)+rho*n; % Received signal
        % ****** 
        % * Rx * 
        % ******
        MSEsol=y./habs;
        [MSEval,MLindex]=min(abs( ones(size(Constellation)) *MSEsol - Constellation*ones(1,SymbolsNum) ).^2,[],1);
        MLindex=MLindex';
        RxBits=BitMapping(MLindex);
        RxBits=int2bit(RxBits',BPS); 
        RxBits=RxBits(:);
        % ********
        % * EVAL *
        % ********
        % Evaluate preformance
        SymErrs(TxNindex,ModulationTypeIndex,RhoIndex)=sum(MLindex~=TxSymbolIndices);
        BitErrs(TxNindex,ModulationTypeIndex,RhoIndex)=sum(RxBits~=TxBits);
        % ******
        % * MC *
        % ******
        KeepGoing=RealizationIndex<MaxRealizations && (SymErrs(TxNindex,ModulationTypeIndex,RhoIndex)<MaxSymbolErrors);
    end
    RealizedSymbolsNum(TxNindex,ModulationTypeIndex,RhoIndex)=RealizationIndex*SymbolsNum;
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
    for TxNindex=1:numel(TxNs)
        LineStyle=num2linestyle(TxNindex);
        Marker=num2marker(TxNindex);
        TxN=TxNs(TxNindex);        
        figure(hSym);        
        semilogy(SNRs, (squeeze(SER(TxNindex,ModulationTypeIndex,:))),strcat(LineStyle, Marker))
        hold on
        figure(hBit)
        semilogy(SNRs, (squeeze(BER(TxNindex,ModulationTypeIndex,:))),strcat(LineStyle,Marker) )
        hold on
    end
    figure(hSym);  
    legend(RXNnames{:})
    title(['MRC: SER vs SNR for modulationType '  ModulationsNames{ModulationType}])
    hold off
    grid on
    xlabel('SNR [dB]')
    ylabel('SER')
    figure(hBit);  
    legend(RXNnames{:})
    title(['MRC: BER vs SNR for modulationType '  ModulationsNames{ModulationType}])
    hold off
    grid on
    xlabel('SNR [dB]')
    ylabel('BER')
end
