function [BER1, BER2] = Tra_NLSE_Rec(sim_type, useEqualizer, samples_to_eyediagram, should_plot)
%%sim_type - if you set this to "signal", the normal PAM signal will be
%transmitted. If you set sim_type to "impulse", Matlab will calculate
%propagation of impulse.

%%useEqualizer - you dont have to use equalizer, but if you want to -
%please set useEqualizer to true. Then BER after equalization will be in
%BER2 var.

%samples_to_eyediagram - in some computers, there may be problems with
%making figure with a lot of data. The Matlab will draw
%samples_to_eyediagram samples in the figure

%should_plot - if you make e.g. for loop to make BER=f(Dispersion) it is
%good to set should_plot to false. Then plotting is turned off.


global nrml2 nrml
global M Baudrate FWHM Wavelength Dispersion DSlope Gamma Length dLength Att RefractiveInex MDF TransmitterCurrent TransmitterBias N
sps=16;%upsampling factor
filter_samples=12;%filter length

Samplerate=Baudrate*sps; %Sample rate
fmax=Samplerate/2;

%Raised cosine filter for preparing signal for transmission in channel.
Rolloff=1;     %roll off factor
raised_cosine_filter=rcosdesign(Rolloff,filter_samples,sps);raised_cosine_filter=raised_cosine_filter(1:end-1); %delete last sample!!

c=3e8;
L_n = Length/dLength; %number of steps in length

wavelength = Wavelength*1e-9;
fwhm = FWHM*1e-12;
FWHM_el = c/wavelength.^2*fwhm;
D = Dispersion*1e-3;
dslope = DSlope*1e6;
beta1 = -RefractiveInex/c/1e3;
beta2 = -D*wavelength^2./(2*pi*c);
beta3 = (wavelength^3/(2*pi*c)^2)*(2*D+wavelength*dslope);
Aeff=pi*(MDF/2)^2;
gamma = Gamma*2*pi/(wavelength*Aeff)/(1e-3);%2.2e-3;%

% MODULATION
if sim_type == "impulse"
    beta1 = 0;
    Baudrate = 1e13;
    sps=1;
    Samplerate = Baudrate*sps;
    N = 2^17;
    impulse_width = 2e-12;
    FWHM_el = 1/impulse_width;
    xx = 1:2^17;xx=xx/Samplerate;
    x_in = exp(-(xx-max(xx)/2).^2/((impulse_width/2)^2));
else
    data=randi(M,1,N)-1;
    sym_send = pammod(data,M,0,'Gray');
    sym_send = sym_send/max(abs(sym_send));
    sym_u=upsample(sym_send,sps);
    x_in=conv(sym_u,raised_cosine_filter);
    x_in=x_in/max(abs(x_in));
end
K = length(x_in);
if mod(K,2)==1
    x_in=[x_in 0];
    K=K+1;
end
f = linspace(0,Samplerate-Samplerate/K,K);
f=[f(1:K/2) -f(K/2+1:-1:2) ];
omega = 2*pi*f;

%%%% TRANSMITER
x_in = TransmitterBias+TransmitterCurrent*x_in;

%phase fluctuations
%10.1109/ISIT.2013.6620632
rnLIM = 1;
phi = zeros(1,length(f(:)'));
transmiter_out = phi;
for i=1:length(f)
    if i>1
        fi = phi(i-1) + randn*2*pi*FWHM_el/Samplerate; %==2*pi*beta*DELTA :)
        phi(i) = fi;
    end
    transmiter_out(i) = exp(1i*phi(i));
end
x_in = sqrt(x_in).*transmiter_out;

%%% CALCULATIONS FOR STEP METHOD from 9780128170434 (Agrawal)
current_L = dLength;
a1_for_full_part = 10^(-Att*current_L/20);
a2_for_full_part = exp(1i*(beta1*current_L*omega+beta2*current_L*(omega.^2)/2+beta3*current_L*omega.^3/6));

current_L = dLength/2;
a1_for_half_part = 10^(-Att*current_L/20);
a2_for_half_part = exp(1i*(beta1*current_L*omega+beta2*current_L*(omega.^2)/2+beta3*current_L*omega.^3/6));

%%Calculate input power
power_IN = rms(x_in)^2;txt=sprintf('%.3f',10*log10(power_IN*1e3));
disp(['Input power: ' txt ' [dBm]'])

%make temp variable
x_signal_TEMP = x_in;

%%%
%<- MAIN NLSE
%%%
if gamma == 0
    x_signal_TEMP_FFT = (fft(x_signal_TEMP));
    for i=1:L_n
        if i ~= 1 && i ~= L_n
            x_signal_TEMP_FFT = ((x_signal_TEMP_FFT.*a2_for_full_part.*a1_for_full_part));
        else
            x_signal_TEMP_FFT = ((x_signal_TEMP_FFT.*a2_for_half_part.*a1_for_half_part));
        end
    end
    x_signal_TEMP = ifft(x_signal_TEMP_FFT);
else
    for i=1:L_n
        x_signal_TEMP_FFT = (fft(x_signal_TEMP));
        if i ~= 1 && i ~= L_n
            current_L = dLength;
            x_signal_TEMP = ifft(x_signal_TEMP_FFT.*a2_for_full_part.*a1_for_full_part);
            x_signal_TEMP = exp(1j*gamma*abs(x_signal_TEMP).^2*current_L).*x_signal_TEMP;
        else
            current_L = dLength/2;
            x_signal_TEMP = ifft(x_signal_TEMP_FFT.*a2_for_half_part.*a1_for_half_part);
            if i==1
                x_signal_TEMP = exp(1j*gamma*abs(x_signal_TEMP).^2*current_L).*x_signal_TEMP;
            end
        end
    end
end
x_signal_TEMP_FFT = fft(x_signal_TEMP);
x_signal_TEMP = ifft(x_signal_TEMP_FFT.*a2_for_full_part.*a1_for_full_part);
%%%
%MAIN NLSE ->
%%%

%%Calculate output power
power_OUT = rms(x_signal_TEMP)^2;txt=sprintf('%.3f',10*log10(power_OUT*1e3));
disp(['Output power: ' txt ' [dBm]'])

%%Simulate photodiode
NEP_linear = 10^(-11);SNR_set = 10*log10(power_OUT/((NEP_linear^2) * fmax));
x_signal_TEMP=x_signal_TEMP+randn(size(x_signal_TEMP))*sqrt(NEP_linear^2*fmax);
x_signal_TEMP=abs(x_signal_TEMP).^2;


%%Receiver
if sim_type == "impulse"
    x_received=x_signal_TEMP;
    figure();plot(xx*1e12, abs(x_in));hold on;yyaxis right;plot(xx*1e12, x_received);
    xlabel('t [ps]')
else
    %
    %FIRST STAGE OF RECEIVING THE SIGNAL
    x_received=x_signal_TEMP-mean(x_signal_TEMP);
    
    %convolute with rcos filter
    x_out=conv(x_received(1:end),raised_cosine_filter);
    %eyediagram(x_out(1:samples_to_eyediagram), sps*2)
    
    %find delay and delete unwanted symbols
    delay2 = finddelay(nrml2(sym_u),nrml2(x_out));
    x_out_delay_cut = x_out(1+delay2:end);
    
    %eyediagram of received signal - check if downsamling starts at good point!
    %eyediagram(x_out_delay_cut(1:samples_to_eyediagram), sps*2)
    
    %downsample data
    received_data = downsample(x_out_delay_cut(1:end),sps);
    % pamdemod(nrml2(received_data)-mean(nrml2(received_data))/2,M)
    
    %omit some symbols at the beginning and at the end of a vector. This helps to prevent some unwanted errors :D
    chosen_symbols = ceil(sps/2):N-sps*4;
    if mod(length(chosen_symbols),2)==1 %make sure that chosen symbols are even
        chosen_symbols = [chosen_symbols chosen_symbols(end)+1];
    end
    
    %%% CALCULATE BER FROM SNR
    evm = comm.EVM;
    %adapt symbols so that they are like sent using ideal transmitter
    symbols_specially_normalized = (TransmitterBias+sym_send*TransmitterCurrent);symbols_specially_normalized=symbols_specially_normalized-mean(symbols_specially_normalized);
    [rmsEVM] = evm(nrml2(symbols_specially_normalized(chosen_symbols))',nrml2(real(received_data(chosen_symbols)))');
    EVM=rmsEVM*0.01;
    SNR1=20*log10(1./EVM);
    BER1=berawgn(SNR1-10*log10(log2(M))-3,'pam',M);
    disp(['BER: ' sprintf('%.3e',BER1)]);
    clear evm
    if should_plot
        figure(101);
        set(gcf,'Position',[300 300 600 400])
        %     subplot(2,1,1)
        %     hold on;
        %     MyEye((sym_u(chosen_symbols(1:samples_to_eyediagram))),sps*2, Samplerate)
        %     title('Ideal signal')
        %     subplot(2,1,2)
        hold on;
        MyEye((x_out(chosen_symbols(100:samples_to_eyediagram+100))),sps*2, Samplerate)
        ylim([-max((x_out(chosen_symbols(100:samples_to_eyediagram+100)))) max((x_out(chosen_symbols(100:samples_to_eyediagram+100))))])
        title('Signal on photodiode')
    end

    
    
    BER = length(find(abs(abs(nrml2(real(received_data(chosen_symbols))))-abs(nrml2(sym_send(chosen_symbols))))>1/M/2))/N;
    
end

if useEqualizer
    cstl=real(pammod(0:1:M-1,M,0,'Gray'));
    if M==2
        eq = lineareq(17,rls(1),cstl); %#ok<CMRLS,CMLRQ>
        eq.RefTap = 1;
        trainlen = 400;
    elseif M==4
        eq = lineareq(17,rls(1),cstl); %#ok<CMRLS,CMLRQ>
        eq.RefTap = 1;
        trainlen = 400;
    elseif M==8
        eq = dfe(17,1,rls(1),cstl); %#ok<CMDFE,CMRLS>
        eq.RefTap = 1;
        trainlen = 200;
        eq = lineareq(13,rls(1),cstl); %#ok<CMRLS,CMLRQ>
        eq.RefTap = 2;
        trainlen = 100;
    end
    toEq = nrml(symbols_specially_normalized(chosen_symbols));
    receivedEq = nrml(real(received_data(chosen_symbols)));
    eqSig = equalize(eq, receivedEq, toEq(1:trainlen)); %#ok<CMEQU>
    latency = finddelay(toEq, eqSig);
    eqSig = eqSig(1+latency+trainlen:end);
    
    N_to_correct = length(eqSig)-ceil(0.25*length(eqSig(:)'));
    % BER after correction
    evm = comm.EVM;
    toSecondCorrectionEq = toEq(1+trainlen:length(toEq)-latency);toSecondCorrectionEq=toSecondCorrectionEq(1:N_to_correct);
    toEvmEq = eqSig(1:N_to_correct);
    delayTo2Evm = finddelay(nrml2(toSecondCorrectionEq(:)), nrml2(toEvmEq(:)));
    [rmsEVM2] = evm(nrml2(toSecondCorrectionEq(1+delayTo2Evm:end))', nrml2(toEvmEq(1:length(toEvmEq)-delayTo2Evm))');
    EVM2=rmsEVM2*0.01;
    SNR2=20*log10(1./EVM2);
    BER2=berawgn(SNR2-10*log10(log2(M))-3,'pam',M);
    disp(['BER after equalizer: ' sprintf('%.3e',BER1)]);
else
    BER2 = 'You have to change useEqualizer to true';
    disp(['BER after equalizer: ' BER2]);
end

end

