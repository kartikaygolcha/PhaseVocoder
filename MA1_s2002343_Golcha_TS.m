%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code Name: Part 1 Matlab Midterm Assessment 
% Q: Basic Time Stretcher
% Developer: Kartikay Golcha
% UUN: s2002343
% Date :10/11/2019
% University: University of Edinburgh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

%read in our WAV file, and store sample rate in Fs
[x,Fs]=audioread('mozart.wav');
x = 0.5*sum(x,2)';                  %stereo to mono
L=length(x);
L_in=L;
%frame length, for now 50ms
time_frame=00.05;
N = round(Fs*time_frame); 
    
H=0.75;                             %Overlap Percentage
HA=round(N-(N*H));                  %actual overlap in bins

%Stretch Coefficient {Q>1 time stretch}, {Q<1 time compress} {Q=1 same length} 
Q=0.5; 


%Defining Window (Hanning)
q=[0:1:N-1];
win=0.5-0.5*cos(2*pi*q/N);

Nf=floor((N+L)/HA);  %Number of frames
%Zero Padding in the start and end  
n1=zeros(1,((Nf-1)*HA-L));
x=[zeros(1,N),x,n1];
L= length(x);                       %Length of total input vector

Hs=floor(Q*HA);                     %Synthesis Hop Length
y=zeros(1,ceil(Q*L));              %output Vector
r=zeros(1,L);                       %for Window interpolation
n=1:1:N;
for i=0:1:Nf-1
    y(1,n+(i*Hs))=y(1,n+(i*Hs))+(2*HA/N)*x(1,n+(i*HA)).*win;
    r(1,n+(i*HA))=r(1,n+(i*HA))+win;
end

%Extracting the final signals 
x=x(1,N+1:L-((Nf-1)*HA-L_in));
y=y(1,N+1:(length(y)-((Nf-1)*HA-L_in)));

%Plotting 

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
    title("Difference in input and output")
end

figure(2);
subplot(2,1,1)
%PLotting the Spectogram INPUT
MA1_s2002343_Golcha_myspec(x',Fs,N,H,"Input");

subplot(2,1,2)
%PLotting the Spectogram OUTPUT
MA1_s2002343_Golcha_myspec(y',Fs,N,H,"Output");

%Playing the output file
%soundsc(y,Fs);