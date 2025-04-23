% Uncoded MRC Symbol Error Rate (SER) & Bit Error Rate (BER) 
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
RxNs=[ 1:4]; % #Rx antennas
%
iif  = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();
RXNnames = arrayfun(@(x) iif(x==1, 'SISO',true ,['MRC 1x' num2str(x) ]),RxNs,'UniformOutput',false);
tol=1e-12;
MaxRealizations=1e1; 
MaxSymbolErrors=MaxRealizations*N; % Max symbols error events number
ModulationsNames={'BPSK','QPSK','8PSK','16 QAM'};
ChannelmodelNames={'AWGN','Rayleigh'};
rhos=10.^(-SNRs/20);
%
[CM,MT,RHO,RXN]=ndgrid(channelModels,ModulationTypes,rhos,RxNs);
ScenarionsNum=numel(RHO);
SymErrs=zeros(numel(RxNs),numel(ModulationTypes),numel(rhos)); % ChannelModelIndex,ModulationTypeIndex
BitErrs=SymErrs;
RealizedSymbolsNum=zeros(size(SymErrs)); % 
for ScenarioIndex=1:ScenarionsNum
    % Progress report
    disp(['Processing scenaio ' num2str(ScenarioIndex) '/' num2str(ScenarionsNum) ' started at ' datestr(now)])
    drawnow
    % Scenario parameters extraction
    rho=RHO(ScenarioIndex);
    [ChannelModelIndex,ModulationTypeIndex,RhoIndex,RxNindex]=ind2sub(size(RHO),ScenarioIndex);
    channelModel=channelModels(ChannelModelIndex);
    ModulationType=ModulationTypes(ModulationTypeIndex);
    RxN=RxNs(RxNindex); % #Rx antennas
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
        s=ones(RxN,1)*Constellation(TxSymbolIndices).';
        % *******
        % * PHY *
        % *******
        n=(randn(RxN,SymbolsNum)+1i*randn(RxN,SymbolsNum))/sqrt(2); % noise
        switch channelModel
            case 1 % AWGN
                h=ones(RxN,SymbolsNum);
            case 2 % Rayleigh
                h=(randn(RxN,SymbolsNum)+1i*randn(RxN,SymbolsNum))/sqrt(2) ; % Flat fading
            otherwise
                error("invalid channel model.")
        end        
        y=h.*s+rho*n; % Received signal
        % ******
        % * Rx * 
        % ******
        h2=ones(RxN,1)*sum(abs(h).^2,1);
        MSEsol=sum(conj(h).*y./h2,1);
        [MSEval,MLindex]=min(abs( ones(size(Constellation)) *MSEsol - Constellation*ones(1,SymbolsNum) ).^2,[],1);
        MLindex=MLindex';
        RxBits=BitMapping(MLindex);
        RxBits=int2bit(RxBits',BPS); 
        RxBits=RxBits(:);
        % ********
        % * EVAL *
        % ********
        % Evaluate preformance
        SymErrs(RxNindex,ModulationTypeIndex,RhoIndex)=SymErrs(RxNindex,ModulationTypeIndex,RhoIndex) + sum(MLindex~=TxSymbolIndices);
        BitErrs(RxNindex,ModulationTypeIndex,RhoIndex)=BitErrs(RxNindex,ModulationTypeIndex,RhoIndex) + sum(RxBits~=TxBits)/BPS;
        % ******
        % * MC *
        % ******
        KeepGoing=RealizationIndex<MaxRealizations && (SymErrs(RxNindex,ModulationTypeIndex,RhoIndex)<MaxSymbolErrors);
    end
    RealizedSymbolsNum(RxNindex,ModulationTypeIndex,RhoIndex)=RealizationIndex*SymbolsNum;
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
    for RxNindex=1:numel(RxNs)
        LineStyle=num2linestyle(RxNindex);
        Marker=num2marker(RxNindex);
        RxN=RxNs(RxNindex);        
        figure(hSym);        
        semilogy(SNRs, (squeeze(SER(RxNindex,ModulationTypeIndex,:))),strcat(LineStyle, Marker))
        hold on
        figure(hBit)
        semilogy(SNRs, (squeeze(BER(RxNindex,ModulationTypeIndex,:))),strcat(LineStyle,Marker) )
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
