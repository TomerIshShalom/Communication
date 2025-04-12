% Uncoded SISO Symbol Error Rate (SER) Monte-Carlo simulation
close all
clear all
clc
drawnow
% Setup params
SNRs=-5:10; % SNRs under consideration [dB]
N=2^12; % Block length
channelModels=[ 1:2]; % 1 - AWGN, 2 - Rayleigh 
ModulationTypes=[2]; % 1 - BPSK, 2 - QPSK, 3- 8PSK, 4 - 16-QAM
%
MinErrorEvents=1e6;
MaxRealizations=1e4;
ModulationsNames={'BPSK','QPSK','8PSK','16 QAM'};
ChannelmodelNames={'AWGN','Rayleigh'};
rhos=10.^(-SNRs/20);
%
[CM,MT,RHO]=ndgrid(channelModels,ModulationTypes,rhos);
ScenarionsNum=numel(RHO);
SymErrs=zeros(numel(channelModels),numel(ModulationTypes),numel(rhos)); % ChannelModelIndex,ModulationTypeIndex
RealizedSymbolsNum=zeros(size(SymErrs)); % 
for ScenarioIndex=1:ScenarionsNum
    rho=RHO(ScenarioIndex);
    [ChannelModelIndex,ModulationTypeIndex,RhoIndex]=ind2sub(size(RHO),ScenarioIndex);
    channelModel=channelModels(ChannelModelIndex);
    ModulationType=ModulationTypes(ModulationTypeIndex);
    %
    switch ModulationType
        case 1 % BPSK
            BPS=1; % Bits Per Symbol
            BitMapping=[1 0]';
            % Constellation design 
            Constellation=(-1).^BitMapping;                
        case 2 % QPSK
            BPS=2; % Bits Per Symbol
            BitMapping=[0 1 3 2]'; % Grey code            
            Constellation=(-1i).^BitMapping;    
        case 3
            BPS=3
            BitMapping=[0 1 3 2 6 7 5 4];
            Constellation=exp(1i*2*pi*(0:7).');
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
    Es=sum(abs(Constellation).^2);
    Constellation=Constellation/sqrt(Es); % Normalize constellation
    % Verify constellation statistics
    Es=sum(abs(Constellation).^2);
    Ms=mean(Constellation);
    disp(['Test constellation normalization of scenario ' num2str(ScenarioIndex) '  - Var=' num2str(Es) ' mean=' num2str(mean(Constellation))])
    %
    RealizationIndex=0;
    KeepGoing=true;
    while KeepGoing
        RealizationIndex=RealizationIndex+1;
        % ******
        % * Tx *
        % ******        
        b=randi(2^BPS,[SymbolsNum,1]);
        s=Constellation(b);
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
        % ********
        % * EVAL *
        % ********
        % Evaluate preformance
        SymErrs(ChannelModelIndex,ModulationTypeIndex,RhoIndex)=sum(MLindex~=b);
        % ******
        % * MC *
        % ******
        KeepGoing=RealizationIndex<=MaxRealizations;
    end
    RealizedSymbolsNum(ChannelModelIndex,ModulationTypeIndex,RhoIndex)=RealizationIndex*SymbolsNum;
end
% Analysis
SER=SymErrs./RealizedSymbolsNum;
for ModulationTypeIndex=1:numel(ModulationTypes)
    ModulationType=ModulationTypes(ModulationTypeIndex );
    figure('Name',ModulationsNames{ModulationType} );
    for ChannelModelIndex=1:numel(channelModels)
        channelModel=channelModels(ChannelModelIndex);
        semilogy(SNRs, (squeeze(SER(ChannelModelIndex,ModulationTypeIndex,:))))
        hold on
    end
    legend(ChannelmodelNames{:})
    title(['ModulationType '  ModulationsNames{ModulationType}])
    hold off
    grid on
    xlabel('SNR [dB]')
    ylabel('SER')
end
