%Fiber transmission
%FiberTeam 2020
%We would like to acknowledge prof. Grzegorz Stepniak. Here, I used a concept
%of Wiener process as random phase fluctuations in transmitter.
clear all;clc;
close all;

global nrml2 nrml
nrml2 = @(x) (x-min(x))/(max(x)-min(x));
nrml = @(x) (x)/max(x);
global M Baudrate FWHM Wavelength Dispersion DSlope Gamma Length dLength Att RefractiveInex MDF TransmitterCurrent TransmitterBias N

M=4; %M-ary PAM e.g. if M = 2, then PAM-2 transmission is used :)
N=6000;%number of transmitted symbols
FWHM = 20; %in picometers 0.08 corresponds to 10 MHz
Baudrate =10e9;
Wavelength = 1550; %in nanometers

Length = 20; %in kilometers
dLength = Length/101;%step of calculations in kilometers!

Att = 0.18; %attenuation coeff of fiber [dB/km] e.g. value: 0.2
Dispersion = 17; %dispersion coeff in ps/(nm*km) e.g. value: 20
DSlope = 0.062; %dispersion slope in ps/(nm^2*km) e.g. value: 0.08
Gamma = 2.6e-20; %Kerr nonlinearity in m^2/W
RefractiveInex = 1.447; %calculated for fiber delay. I encourage you to set it to 0, because we dont send appropriate number of symbols to find the delay in the data!
MDF = 10.36e-6; %mode field diameter e.g. value: 10.4e-6

TransmitterCurrent = 10e-3;
TransmitterBias = 10e-3;

type = "signal";
%type = "impulse";
samples_to_eyediagram = 4000;
use_equalizer = false;

%PLEASE, read first lines of comments in Tra_NLSE_Rec!!
Tra_NLSE_Rec(type, use_equalizer, samples_to_eyediagram, true);

% %
% BER_s = [];Dispersion_s=[20:1:90];
% for Dispersion=Dispersion_s
%     [B1,~]=Tra_NLSE_Rec(type, use_equalizer, samples_to_eyediagram,false);
%     BER_s = [BER_s, B1];
% end
% figure();plot(Dispersion_s, log10(BER_s), 'linewidth',3);ylabel('log10(BER)');xlabel('Dispersion [ps/(nm*km)]');ylim([-20 0]);grid on;

%
BER_s = [];
M_s = [2 4 8];
Length_s=[10:5:100];
itL = 1;
for Length=Length_s
    for i=1:length(M_s)
        M = M_s(i);
        %Length = 20; %in kilometers
        dLength = Length/101;%step of calculations in kilometers!
        [B1,~]=Tra_NLSE_Rec(type, use_equalizer, samples_to_eyediagram,false);
        BER_s(itL,i) = B1;
    end
    itL = itL + 1;
end

figure('color','w');plot(Length_s, log10(BER_s), 'o', 'linewidth',3);ylabel('log10(BER)');xlabel('Długość [km]');ylim([-10 0]);grid on;
legend({'PAM2','PAM4','PAM8'},'box','off','Location','best')