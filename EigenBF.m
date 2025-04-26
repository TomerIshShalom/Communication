% Uncoded Eigen beamforming (with one Rx antenna)
% Symbol Error Rate (SER) & Bit Error Rate (BER) 
% Monte-Carlo simulation
close all
clear all
clc
drawnow
% Setup params
SNRs=0:30; % SNRs under consideration [dB]
N=2*50000; % Block length
channelModels=[ 2]; % 1 - AWGN, 2 - Rayleigh 
ModulationTypes=[2]; % 1 - BPSK, 2 - QPSK, 3- 8PSK, 4 - 16-QAM
TxNs=[ 2]; % #Rx antennas
RxNs=[ 2,4]; % #Rx antennas
%
tol=1e-12;
MaxRealizations=1e0; 
MaxSymbolErrors=MaxRealizations*N; % Max symbols error events number
ModulationsNames={'BPSK','QPSK','8PSK','16 QAM'};
ChannelmodelNames={'AWGN','Rayleigh'};
rhos=10.^(-SNRs/20);
%
[CM,MT,RHO,TXN,RXN]=ndgrid(channelModels,ModulationTypes,rhos,TxNs,RxNs);
ScenarionsNum=numel(RHO);
SymErrs=zeros(numel(TxNs),numel(RxNs),numel(ModulationTypes),numel(rhos)); % ChannelModelIndex,ModulationTypeIndex
BitErrs=SymErrs;
RealizedSymbolsNum=zeros(size(SymErrs)); % 
for ScenarioIndex=1:ScenarionsNum
    % Progress report
    disp(['Processing scenaio ' num2str(ScenarioIndex) '/' num2str(ScenarionsNum) ' started at ' datestr(now)])
    drawnow
    % Scenario parameters extraction
    rho=RHO(ScenarioIndex);
    [ChannelModelIndex,ModulationTypeIndex,RhoIndex,TxNindex,RxNindex]=ind2sub(size(RHO),ScenarioIndex);
    channelModel=channelModels(ChannelModelIndex);
    ModulationType=ModulationTypes(ModulationTypeIndex);
    TxN=TxNs(TxNindex); % #Tx antennas
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
    MLindex=zeros(SymbolsNum,1);
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
        s=Constellation(TxSymbolIndices).';
        % Feedback 
        switch channelModel
            case 1 % AWGN
                h=ones(RxN,TxN,SymbolsNum);
            case 2 % Rayleigh
                h=(randn(RxN,TxN,SymbolsNum)+1i*randn(RxN,TxN,SymbolsNum))/sqrt(2) ; % Flat fading
            otherwise
                error("invalid channel model.")
        end   
        for SymbolIndex=1:SymbolsNum
            % Precoding  
            H=squeeze(h(:,:,SymbolIndex));
            [U,S,V]=svd(H);
            w=V(:,1); % Singular vector            
            MRCvec=(H*w)'/S(1,1)^2;
            ps=w*s(SymbolIndex);
            % *******
            % * PHY *
            % *******
            n=(randn(RxN,1)+1i*randn(RxN,1))/sqrt(2); % noise    
            y=H*ps+rho*n; % Received signal
            % ****** 
            % * Rx * 
            % ******
            MSEsol=MRCvec*y;
            [MSEval,MLindex(SymbolIndex)]=min(abs( ones(size(Constellation)) *MSEsol - Constellation ).^2,[],1);            
        end
        RxBits=BitMapping(MLindex);
        RxBits=int2bit(RxBits',BPS); 
        RxBits=RxBits(:);
        % ********
        % * EVAL *
        % ********
        % Evaluate preformance
        SymErrs(TxNindex,RxNindex,ModulationTypeIndex,RhoIndex)=SymErrs(TxNindex,RxNindex,ModulationTypeIndex,RhoIndex)+sum(MLindex~=TxSymbolIndices);
        BitErrs(TxNindex,RxNindex,ModulationTypeIndex,RhoIndex)=BitErrs(TxNindex,RxNindex,ModulationTypeIndex,RhoIndex)+sum(RxBits~=TxBits)/BPS;
        % ******
        % * MC *
        % ******
        KeepGoing=RealizationIndex<MaxRealizations && (SymErrs(TxNindex,ModulationTypeIndex,RhoIndex)<MaxSymbolErrors);
    end
    RealizedSymbolsNum(TxNindex,RxNindex,ModulationTypeIndex,RhoIndex)=RealizationIndex*SymbolsNum;
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
    [TXN,RXN]=ndgrid(TxNs,RxNs);
    ScenarionsNum=numel(TXN);
    legends={};
    for ScenarioIndex=1:ScenarionsNum        
        [TxNindex,RxNindex]=ind2sub(size(TXN),ScenarioIndex);
        TxN=TxNs(TxNindex); % #Tx antennas
        RxN=RxNs(RxNindex); % #Rx antennas
        LineStyle=num2linestyle(TxNindex);
        Marker=num2marker(RxNindex);        
        %
        figure(hSym);        
        semilogy(SNRs, (squeeze(SER(TxNindex,RxNindex,ModulationTypeIndex,:))),strcat(LineStyle, Marker))
        hold on
        figure(hBit)
        semilogy(SNRs, (squeeze(BER(TxNindex,RxNindex,ModulationTypeIndex,:))),strcat(LineStyle,Marker) )
        hold on
        legends{ScenarioIndex}=[num2str(TxN) 'x' num2str(RxN)];
    end
    figure(hSym);  
    legend(legends{:})
    title(['MRC: SER vs SNR for modulationType '  ModulationsNames{ModulationType}])
    hold off
    grid on
    xlabel('SNR [dB]')
    ylabel('SER')
    figure(hBit);  
    legend(legends{:})
    title(['MRC: BER vs SNR for modulationType '  ModulationsNames{ModulationType}])
    hold off
    grid on
    xlabel('SNR [dB]')
    ylabel('BER')
end
