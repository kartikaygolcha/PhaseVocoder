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
[x,Fs]=audioread('mozart.wav');
x = 0.5*sum(x,2)';                  %stereo to mono
x2=x;
L_in=length(x);

%frame length, for now 50ms
time_frame=0.05;
N = round(Fs*time_frame); 
fft_size=ceil(log(N)/log(2));
    
H=0.75;                             %Overlap Percentage
HA=round(N-(N*H));                  %actual overlap in bins

%Stretch Coefficient {Q>1 time stretch}, {Q<1 time compress} {Q=1 same length} 
Q=1; 

%Defining Window (Hanning)
q=[0:1:N-1];
win=0.5-0.5*cos(2*pi*q/N);
Nf=floor((N+L_in)/HA);
Nfft=2^fft_size;

%Zero Padding in the start and end  
if(((Nf-1)*HA-L_in)>N/2)
    n1=(Nf-1)*HA-L_in+N;
    n2=N;
else 
    n1=(Nf-1)*HA-L_in;
    n2=0;
end    
x=[zeros(1,N),x,zeros(1,n1)];
L= length(x);                       %reinitialize Length of total input vector
 
Hs=floor(Q*HA);                     %Synthesis Hop Length
y=zeros(1,ceil(n2+N+(Hs*(Nf-1))));               %output Vector
r=zeros(1,L);                       %for Window interpolation
n=1:1:N;

%frequency vector for N-length DFT
phi_m=zeros(1,2^fft_size);
theta=zeros(1,2^fft_size);
wmk=zeros(1,2^fft_size);
bin_freq=2*pi*linspace(1,2^fft_size,2^fft_size)/(2^fft_size);

for i=0:1:Nf
    fft_in=fft((x(1,n+i*HA)),2^fft_size);
    Xmag=abs(fft_in);
    Xang=angle(fft_in);
    wmk=ppa((Xang-phi_m-(bin_freq.*HA))/HA);
    theta=theta+(wmk*Hs)+(bin_freq.*Hs);
    Y=Xmag.*exp(j*theta);
    
    %Check Hermitian Symmetry
    Err=Y(2:Nfft/2)-(flip(conj(Y((2^fft_size/2 + 2):Nfft))));
    assert(abs(real(Y(1)))<=abs(complex(Y(1))) || Y(1)~=0);
    %assert(abs(real(Y(2049)))<=abs(complex(Y(2049))) || Y(2049)~=0);
    assert(sum(Err)<0.01);
    
    x_ifft=ifft(Y,'symmetric');
    x_ifft=x_ifft(1,1:N).*win;
    y(1,n+(i*Hs))=y(1,n+(i*Hs))+((2*HA/N)*x_ifft);
    r(1,n+(i*HA))=r(1,n+(i*HA))+win;
    phi_m=Xang;
end
%Extracting the final signals 
x=x(1,N+1:L-n1);
y=y(1,N+1:(length(y)-((Nf-1)*HA-L_in))-n2);

%Plotting 

%soundsc(y,Fs);
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
    plot(r);
    xlabel("Time");
    ylabel("Amplitude");
    title("Window Interpolation");
    
    subplot(4,1,4); 
    plot(x-y);
    xlabel("Time");
    ylabel("Amplitude");
    title("Difference in input and output");
end
figure(2);
subplot(2,1,1)
%PLotting the Spectogram INPUT
MA1_s2002343_Golcha_myspec(x',Fs,N,H,"Spectogram for input");

subplot(2,1,2)
%PLotting the Spectogram OUTPUT
MA1_s2002343_Golcha_myspec(y',Fs,N,H,"Spectogram for output");

%Playing the output file
%soundsc(y,Fs);