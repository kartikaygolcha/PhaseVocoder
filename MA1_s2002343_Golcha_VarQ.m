%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code Name: Part 1 Matlab Midterm Assessment 
% Q: Phase Vocoder Time Stretcher
% Developer: Kartikay Golcha
% UUN: s2002343
% Date :10/11/2019
% University: University of Edinburgh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
ppa=inline('mod(a+pi,2*pi)-pi','a'); %Function to WrapToPi

%read in our WAV file, and store sample rate in Fs
[x,Fs]=audioread('piano.wav');
x = 0.5*sum(x,2)';                  %stereo to mono

%frame length, for now 50ms
time_frame=0.05;
N = round(Fs*time_frame); 
fft_size=ceil(log(N)/log(2));
    
H=0.75;                             %Overlap Percentage
HA=round(N-(N*H));                  %actual overlap in bins

%Stretch Coefficient {Q>1 time stretch}, {Q<1 time compress} {Q=1 same length} 


%Defining Window (Hanning)
q=[0:1:N-1];
win=0.5-0.5*cos(2*pi*q/N);

%Zero Padding in the start and end  
x=[zeros(1,N),x,zeros(1,round(N))];
L= length(x);                       %Length of total input vector

Nf=floor((L-N)/HA);                 %Number of frames
Q=exp(linspace(1.2,1/2,Nf))/2.73; 
Hs=floor(Q*HA);                     %Synthesis Hop Length
hh=sum(Q.*HA));
y=zeros(1,ceil(hh));              %output Vector
r=zeros(1,L);                       %for Window interpolation
n=1:1:N;


%frequency vector for N-length DFT
phi_m=zeros(1,2^fft_size);
theta=zeros(1,2^fft_size);
wmk=zeros(1,2^fft_size);
bin_freq=2*pi*linspace(1,2^fft_size,2^fft_size)/(2^fft_size);

for i=1:1:Nf
    fft_in=fft((x(1,n+i*HA)),2^fft_size);
    Xmag=abs(fft_in);
    Xang=angle(fft_in);
    wmk=ppa((Xang-phi_m-(bin_freq.*HA))/HA);
    theta=theta+(wmk*Hs(i))+(bin_freq.*Hs(i));
    Y=Xmag.*exp(j*theta);
    x_ifft=ifft(Y,'symmetric');
    x_ifft=x_ifft(1,1:N).*win;
    y(1,n+(i*Hs(i)))=y(1,n+(i*Hs(i)))+((1/2)*x_ifft);
    r(1,n+(i*HA))=r(1,n+(i*HA))+win;
    phi_m=Xang;
end

%Extracting the final signals 
x=x(1,N:L-round(2*N-HA));
y=y(1,N:(length(y)-round(2*N-HA)));

%Plotting 
soundsc(y,Fs);
subplot(4,1,1);
plot(x);
xlabel("Time");
ylabel("Amplitude");
title("Input of Time Stretcher");

subplot(4,1,2);
plot(y);
xlabel("Time");
ylabel("Amplitude");
title("Output of Time Stretcher");
if Q==1
    subplot(4,1,3);
    MA1_s2002343_Golcha_myspec(y',Fs,N,H);
    xlabel("Time");
    ylabel("Amplitude");
    title("Window Interpolation");
    
    subplot(4,1,4); 
    plot(x-y);
    xlabel("Time");
    ylabel("Amplitude");
    title("Difference in input and output")
end

%Playing the output file
%soundsc(y,Fs);