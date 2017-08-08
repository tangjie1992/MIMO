clc
clear
close all

N_frame=100;  N_packet=100;
b=2; M=2^b;  % Number of bits per symbol and Modulation order
mod_obj = modem.qammod('M',M,'SymbolOrder','gray','InputType','bit');
demod_obj = modem.qamdemod(mod_obj);
NT=4;  NR=4;  sq2=sqrt(2);  I=eye(NR); I2=eye(2*NR);
N_pbits = N_frame*NT*b;
N_tbits = N_pbits*N_packet;

SNRdBs = [0:4:40];
for i_SNR=1:length(SNRdBs)
   SNRdB = SNRdBs(i_SNR); 
   noise_var = NT*0.5*10^(-SNRdB/10); 
   sigma = sqrt(noise_var);
   rand('seed',1); randn('seed',1);  
   N_ebits = 0; N_ebits2 = 0;
   %%%%%%%%%%%%% Transmitter %%%%%%%%%%%%%%%%%%
   for i_packet=1:N_packet
      H = (randn(NR,NT)+1i*randn(NR,NT))/sq2;
      msg_bit = randint(N_pbits,1); % bit generation
      symbol = modulate(mod_obj,msg_bit).';
      Scale = modnorm(symbol,'avpow',1); % normalization
      Symbol_nomalized = reshape(Scale*symbol,NT,N_frame); 
      
      Tx_signal = Symbol_nomalized;   
      %%%%%%%%%%%%% Channel and Noise %%%%%%%%%%%%% 
      noise = sigma*(randn(NR,N_frame)+1i*randn(NR,N_frame));
      Rx_signal = H*Tx_signal+noise;
      %%%%%%%%%%%%%% Receiver %%%%%%%%%%%%%%%%%%%%%
      % method 1
      Symbol_hat = inv(H'*H+noise_var*I)*H'*Rx_signal;
      msg_hat = demodulate(demod_obj,Symbol_hat(:)/Scale);
      N_ebits = N_ebits + sum(msg_hat~=msg_bit);
      % method 2
      x_hat = LRAD_MMSE2(H,Rx_signal,noise_var,M,0.75);
      msg_hat2 = demodulate(demod_obj,x_hat(:)/Scale);
      N_ebits2 = N_ebits2 + sum(msg_hat2~=msg_bit);
   end
   BER(i_SNR,1) = N_ebits/N_tbits;
   BER(i_SNR,2) = N_ebits2/N_tbits;
end
figure;
semilogy(SNRdBs,BER(:,1),'-.k^'); hold on; 
semilogy(SNRdBs,BER(:,2),'-r+');  hold off;
grid on;xlabel('SNR[dB]'), ylabel('BER');
title('MIMO Detection');legend('MMSE','MMS+LLL')
