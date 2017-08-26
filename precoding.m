clc
clear
close all

N_frame=100;  N_packet=100;
b=2; M=2^b;  % Number of bits per symbol and Modulation order
mod_obj = modem.qammod('M',M,'SymbolOrder','gray','InputType','bit');
demod_obj = modem.qamdemod(mod_obj);
NT=4;  NR=4;  sq2=sqrt(2);  I=eye(NR); I2=eye(2*NR);
N_pbits = N_frame*NT*b; N_tbits = N_pbits*N_packet;
SNRdBs = [30:4:50];
for i_SNR=1:length(SNRdBs)
   SNRdB = SNRdBs(i_SNR); 
   noise_var = NT*0.5*10^(-SNRdB/10); 
   sigma = sqrt(noise_var);
   rand('seed',1); randn('seed',1);  
   N_ebits = 0; N_ebits2 = 0;
   
   for i_packet=1:N_packet
      H = (randn(NR,NT)+1i*randn(NR,NT))/sq2;
      msg_bit = randint(N_pbits,1); % bit generation
      symbol = modulate(mod_obj,msg_bit).';
      Scale = modnorm(symbol,'avpow',1); % normalization
      Symbol_nomalized = reshape(Scale*symbol,NT,N_frame); 
      %%%%%%%%%%%%% Transmitter %%%%%%%%%%%%%%%%%%
      W_zf = H'*inv(H*H');
      beta_zf = sqrt(NT/trace(W_zf*W_zf'));
      Tx_zf = beta_zf*W_zf*Symbol_nomalized;
      
      W_mmse = H'*inv(H*H'+noise_var*I);
      beta_mmse = sqrt(NT/trace(W_mmse*W_mmse'));
      Tx_mmse = beta_mmse*W_mmse*Symbol_nomalized;
      %%%%%%%%%%%%% Channel and Noise %%%%%%%%%%%%% 
      noise = sigma*(randn(NR,N_frame)+1i*randn(NR,N_frame));
      Rx_zf = H*Tx_zf+noise;
      Rx_mmse = H*Tx_mmse+noise;
      %%%%%%%%%%%%%% Receiver %%%%%%%%%%%%%%%%%%%%%
      Symbol_hat_zf = Rx_zf/beta_zf;
      msg_hat_zf = demodulate(demod_obj,Symbol_hat_zf(:)/Scale);
      N_ebits = N_ebits + sum(msg_hat_zf ~= msg_bit);
      
      Symbol_hat_mmse = Rx_mmse/beta_mmse;
      msg_hat_mmse = demodulate(demod_obj,Symbol_hat_mmse(:)/Scale);
      N_ebits2 = N_ebits2 + sum(msg_hat_mmse ~= msg_bit);
   end
   BER(i_SNR,1) = N_ebits/N_tbits;
   BER(i_SNR,2) = N_ebits2/N_tbits;
end
figure;semilogy(SNRdBs,BER(:,1),'-*');hold on
semilogy(SNRdBs,BER(:,2),'-o');hold off
legend('ZF-precoding','MMSE-precoding')
grid on;
% test github
