% Replicating
% https://www.usna.edu/NAOE/_files/documents/Courses/EN455/AY20_Notes/EN455LaboratoryAssignmentsAY20_IrregularWaves.pdf 

% Linear supersition - adding together regular waves to model irregular
% waves
close all
clear all



Fs = 100; % sampling frequency
SuperWave1 = zeros(1,Fs);
for i = 1:4 % number of regular waves
    t = linspace(0,10,Fs);
    freqHz = rand*2;
    freq = 2*pi*freqHz;
    
    freqHz2 = rand*2;
    freq2 = 2*pi*freqHz2;
    amp1 = rand;
    amp2 = rand;
    Wave1 = amp1*cos(freq*t);
    Wave2 = amp2*cos(freq2*t+pi/4);%pi/4 is a phase shift

    SuperWave1 = SuperWave1 + Wave1 + Wave2;
end

%conFigure(11,3)
f1 = figure;
plot(t,SuperWave1)
xlabel('Time [s]')
ylabel('Wave elevation [m]')




% Perform Fourier transform
IrregularWave1 = SuperWave1;
L = length(IrregularWave1);

if mod(L,2) == 0
    n = L/2;
else
    n = (L-1)/2;
end

WaveFFT = fft(IrregularWave1);
WaveFFt = WaveFFT/L;

WaveMag2 = abs(WaveFFT);
WavePhase2 = angle(WaveFFT)*180/pi;
WaveMag = WaveMag2(1:n+1);
WaveMag(2:end-1) = 2*WaveMag(2:end-1);
WavePhase = WavePhase2(1:n+1);

f = Fs*(0:n)/L;

f2 = figure;
plot(f,WaveMag)
f3 = figure;
plot(f,WavePhase)
