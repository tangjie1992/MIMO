r = [0:0.01:10];
EbN0 = 10*log10((2.^r-1)./r);

figure;plot(EbN0,r,'-ko');
xlabel('E_b/N_0(dB)');
ylabel('bit/(s*Hz)');
title('Channel Capacity');
grid on
