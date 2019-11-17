%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code Name: Matlab Assessment 1 
% Q: Spectogram
% Developer: Kartikay Golcha
% UUN: s2002343
% Date :10/11/2019
% University: University of Edinburgh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function my_spectogram_out=MA1_s2002343_Golcha_myspec_1(ui,Input,Fs,Sample_window,Overlap)
    x=Input;
    N=Sample_window;
    H=Overlap;
    HA=round(N-(N*H));                  %actual overlap in bins 
    fft_size=ceil(log(N)/log(2));
  
    %Defining Window (Hanning)
    q=[0:1:N-1];
    win=0.5-0.5*cos(2*pi*q/N);

    %Zero Padding in the start and end  
    x=[zeros(N,1);x;zeros(round(N-HA),1)];
    L= length(x);                       %Length of total input vector

    Nf=floor((L-N)/HA);                 %Number of frames
    n=1:1:N;
    fft_in=zeros(2^fft_size,1);
    yfft=zeros(2^(fft_size-1),1);
    f=(Fs./(2*(1:2^(fft_size))));
    t=[1:1:Nf+1];
    for i=1:1:Nf
        fft_in=fft(x(n+(i-1)*HA,1).*win',2^fft_size);
        Xmag=abs(fft_in);
        Xmag=Xmag(1:length(fft_in)/2,1);
        Xang=angle(fft_in);
        yfft=[yfft,Xmag];
    end
    yfft=flip(20*log10(yfft));
    
    f=f(1:2^(fft_size-1));
    imagesc(ui,'XData',t,'YData',f','CData',yfft);
    %imagesc(ui,t,f',yfft);
    colorbar(ui);
    caxis(ui,[-60 0]);
end