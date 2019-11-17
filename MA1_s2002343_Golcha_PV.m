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
if(((Nf-1)*HA-L_in)>N/2)  % zero padding at the end to accouont of the slope of window
    n1=(Nf-1)*HA-L_in+N; 
    n2=N;
else 
    n1=(Nf-1)*HA-L_in;
    n2=0;
end    

x=[zeros(1,N),x,zeros(1,n1)];
L= length(x);                       %reinitialize Length of total input vector
 
Hs=floor(Q*HA);                     %Synthesis Hop Length
y=zeros(1,ceil(n2+N+(Hs*(Nf-1))));  %output Vector
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
    assert(abs(real(Y(Nfft/2+1)))<=abs(complex(Y(Nfft/2 + 1))) || Y(Nfft/2 + 1)~=0);
    assert(sum(Err)<0.01);
    
    x_ifft=ifft(Y,'symmetric');
    x_ifft=x_ifft(1,1:N).*win;
    y(1,n+(i*Hs))=y(1,n+(i*Hs))+((2*HA/N)*x_ifft);
    r(1,n+(i*HA))=r(1,n+(i*HA))+win;
    phi_m=Xang;
end
%Extracting the final signals 
x_in=x(1,N+1:L-n1);
y_out=y(1,N+1:(length(y)-((Nf-1)*HA-L_in))-n2);

t_x=[0:1:length(x_in)-1]/Fs;
t_y=[0:1:length(y_out)-1]/Fs;
t_r=[0:1:length(r)-1]/Fs;
%Plotting 
fig1=figure('Name',"Plots");
subplot(4,1,1);
plot(t_x,x_in);
xlabel("Time(s)");
ylabel("Amplitude");
title("Input of Time Stretcher");

subplot(4,1,2);
plot(t_y,y_out);
xlabel("Time(s)");
ylabel("Amplitude");
title("Output of Time Stretcher");
if Q==1
    subplot(4,1,3);
    plot(t_r,r);
    xlabel("Time(s)");
    ylabel("Amplitude");
    title("Window Interpolation");
    
    subplot(4,1,4); 
    plot(t_x,x_in-y_out);
    xlabel("Time(s)");
    ylabel("Amplitude");
    title("Difference in input and output");
end

%Spectogram Plots
figure('Name',"Spectograms");
subplot(2,1,1)
%PLotting the Spectogram INPUT
MA1_s2002343_Golcha_myspec(x',Fs,N,H,"Spectogram for input");

subplot(2,1,2)
%PLotting the Spectogram OUTPUT
MA1_s2002343_Golcha_myspec(y',Fs,N,H,"Spectogram for output");

%Playing the output file
%soundsc(y,Fs);