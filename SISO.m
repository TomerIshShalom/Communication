% Uncoded SISO Symbol Error Rate (SER) Monte-Carlo simulation
close all
%clear all
%clc
drawnow
% Setup params
SNRs=0:30; % SNRs under consideration [dB]
N=2*500000;%2^12; % Block length
channelModels=[ 1:2]; % 1 - AWGN, 2 - Rayleigh 
ModulationTypes=[2]; % 1 - BPSK, 2 - QPSK, 3- 8PSK, 4 - 16-QAM
%
tol=1e-12;
MaxRealizations=1e1;
MaxSymbolErrorEvents=MaxRealizations*N;
ModulationsNames={'BPSK','QPSK','8PSK','16 QAM'};
ChannelmodelNames={'AWGN','Rayleigh'};
rhos=10.^(-SNRs/20);
%
[CM,MT,RHO]=ndgrid(channelModels,ModulationTypes,rhos);
ScenarionsNum=numel(RHO);
SymErrs=zeros(numel(channelModels),numel(ModulationTypes),numel(rhos)); % ChannelModelIndex,ModulationTypeIndex
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
            % Constellation=(-1i).^BitMapping;  
            Constellation=exp(-1i*(1.5:-1:-1.5)*2*pi/4).';
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
    if abs(Ms)>tol || abs(sqrt(Es)-1)>tol
        Constellation=(Constellation-Ms)/sqrt(Es); % Normalize constellation
        % Verify constellation statistics
        Es=mean(abs(Constellation).^2);
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
        s=Constellation(TxSymbolIndices);
        % *******
        % * PHY *
        % *******
        n=(randn(SymbolsNum,1)+1i*randn(SymbolsNum,1))/sqrt(2); % noise
        switch channelModel
            case 1 % AWGN
                h=ones(SymbolsNum,1);
            case 2 % Rayleigh
                h=(randn(SymbolsNum,1)+1i*randn(SymbolsNum,1))/sqrt(2) ; % Flat fading
            otherwise
                error("invalid channel model.")
        end        
        y=h.*s+rho*n; % Received signal
        % ******
        % * Rx * 
        % ******
        [MSEval,MLindex]=min(abs( (y./h)*ones(size(Constellation'))  - ones(size(y))*(Constellation.') ).^2,[],2);
        RxBits=BitMapping(MLindex);
        RxBits=int2bit(RxBits',BPS); 
        RxBits=RxBits(:);
        % ********
        % * EVAL *
        % ********
        % Evaluate preformance
        SymErrs(ChannelModelIndex,ModulationTypeIndex,RhoIndex)=SymErrs(ChannelModelIndex,ModulationTypeIndex,RhoIndex)+sum(MLindex~=TxSymbolIndices);
        BitErrs(ChannelModelIndex,ModulationTypeIndex,RhoIndex)=BitErrs(ChannelModelIndex,ModulationTypeIndex,RhoIndex)+sum(RxBits~=TxBits)/BPS;
        % ******
        % * MC *
        % ******
        KeepGoing=(RealizationIndex<MaxRealizations) && (SymErrs(ChannelModelIndex,ModulationTypeIndex,RhoIndex)<MaxSymbolErrorEvents);
    end
    RealizedSymbolsNum(ChannelModelIndex,ModulationTypeIndex,RhoIndex)=RealizationIndex*SymbolsNum;
end
% Analysis
SER=SymErrs./RealizedSymbolsNum;
BER=BitErrs./RealizedSymbolsNum;

for ModulationTypeIndex=1:numel(ModulationTypes)
    ModulationType=ModulationTypes(ModulationTypeIndex );
    hSym=figure('Name',ModulationsNames{ModulationType} );
    hBit=figure('Name',ModulationsNames{ModulationType} );
    for ChannelModelIndex=1:numel(channelModels)
        channelModel=channelModels(ChannelModelIndex);
        figure(hSym);        
        semilogy(SNRs, (squeeze(SER(ChannelModelIndex,ModulationTypeIndex,:))))
        hold on
        figure(hBit)
        semilogy(SNRs, (squeeze(BER(ChannelModelIndex,ModulationTypeIndex,:))))
        hold on
    end
    figure(hSym);  
    legend(ChannelmodelNames{:})
    title(['SER vs SNR for modulationType '  ModulationsNames{ModulationType}])
    hold off
    grid on
    xlabel('SNR [dB]')
    ylabel('SER')
    figure(hBit);  
    legend(ChannelmodelNames{:})
    title(['BER vs SNR for modulationType '  ModulationsNames{ModulationType}])
    hold off
    grid on
    xlabel('SNR [dB]')
    ylabel('BER')
end
