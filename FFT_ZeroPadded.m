function [xdft,freq, p ] = FFT_ZeroPadded(In_Signal,PaddingPts, Fs, flag)

lpad = PaddingPts*length(In_Signal);
xdft = fft(In_Signal,lpad);
xdft = xdft(1:lpad/2+1);
xdft = xdft/length(In_Signal);
xdft(2:end-1) = 2*xdft(2:end-1);
freq = 0:Fs/lpad:Fs/2;


if flag == 1
    
   p = plot(freq, abs(xdft),'b');
   title('FFT')
   xlabel('Freq')
   ylabel('Magnitude')
   grid on
  xlim([0 10])

else
    p = NaN;
    
end

end

