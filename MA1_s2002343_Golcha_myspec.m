%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code Name: Matlab Assessment 1 
% Q: Spectogram
% Developer: Kartikay Golcha
% UUN: s2002343
% Date :10/11/2019
% University: University of Edinburgh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function my_spectogram_out=MA1_s2002343_Golcha_myspec(Input,Fs,Sample_window,Overlap,Graph)
    x=Input;
    N=Sample_window;
    H=Overlap;
    HA=round(N-(N*H));                  %actual overlap in bins 
    fft_size=ceil(log(N)/log(2));
     L= length(x); 
  
    %Defining Window (Hanning)
    q=[0:1:N-1];
    win=0.5-0.5*cos(2*pi*q/N);
    Nf=floor((L+N)/HA);                 %Number of frames
   
    %Zero Padding in the start and end  
    x=[zeros(N,1);x;zeros((Nf-1)*HA-L,1)];
    L= length(x);                       %Length of total input vector
       
    n=1:1:N;
    fft_in=zeros(2^fft_size,1);
    yfft=zeros(2^(fft_size-1),1);
    f=((1:2^fft_size)./2^fft_size)*Fs;
    t=[1:1:Nf+1]*HA/Fs;
    for i=1:1:Nf
        fft_in=fft(x(n+(i-1)*HA,1).*win',2^fft_size);
        Xmag=abs(fft_in);
        Xmag=Xmag(1:length(fft_in)/2,1);
        Xang=angle(fft_in);
        Xmag=(20*log10(Xmag));
        yfft=[yfft,Xmag];
        
    end
    f=f(1:2^(fft_size-1));
    title_spec=strcat(Graph," N=",num2str(N/Fs)," Overlap=",num2str(Overlap*100)); 
    xlabel("Time");
    ylabel("Frequency");
    imagesc('XData',t,'YData',f','CData',yfft);
    title(title_spec);
    fft_max=max(yfft(:));
    colorbar;
    caxis([fft_max-60 fft_max]);
end