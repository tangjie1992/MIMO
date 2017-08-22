clc
clear
close all

N_frame=150;  N_packet=10;
b=2; M=2^b;  % Number of bits per symbol and Modulation order
mod_obj = modem.qammod('M',M,'SymbolOrder','gray','InputType','bit');
demod_obj = modem.qamdemod(mod_obj);
NT=4;  NR=4;  sq2=sqrt(2);  I=eye(NR); I2=eye(2*NR);
N_pbits = N_frame*NT*b; N_tbits = N_pbits*N_packet;

load matt
tao = 2*sqrt(2);

SNRdBs = [0:4:40];
for ii=1:length(SNRdBs)
    SNRdB = SNRdBs(ii)
    noise_var = NT*0.5*10^(-SNRdB/10); 
    sigma = sqrt(noise_var);
    rand('seed',1); randn('seed',1);  
    N_ebits_zf = 0; N_ebits_vp = 0;
    for i_packet=1:N_packet
       H = (randn(NR,NT)+1i*randn(NR,NT))/sq2;
       msg_bit = randint(N_pbits,1); % bit generation
       symbol = modulate(mod_obj,msg_bit).';
       Scale = modnorm(symbol,'avpow',1); % normalization
       Symbol_nomalized = reshape(Scale*symbol,NT,N_frame); 
       %%%%%%%%%%%%% Transmitter %%%%%%%%%%%%%%%%%%
       [Tx_signal_zf,beta_zf]=precoding_zf(H,Symbol_nomalized);
       [Tx_signal_vp,beta_vp,arr2]=precoding_vp(H,Symbol_nomalized,matt);
       %%%%%%%%%%%%% Channel and Noise %%%%%%%%%%%%% 
       noise = sigma*(randn(NR,N_frame)+1i*randn(NR,N_frame));
       Rx_signal_zf = H*Tx_signal_zf+noise;
       Rx_signal_vp = H*Tx_signal_vp+noise;
       %%%%%%%%%%%%%% Receiver %%%%%%%%%%%%%%%%%%%%%
       Symbol_hat_zf = Rx_signal_zf/beta_zf;
       msg_hat_zf = demodulate(demod_obj,Symbol_hat_zf(:)/Scale);
       N_ebits_zf = N_ebits_zf + sum(msg_hat_zf ~= msg_bit);
       
       Rx_signal_vp2 = Rx_signal_vp./repmat(beta_vp,NT,1);
       re = mod(real(Rx_signal_vp2)+tao/2,tao)-tao/2;
       im = mod(imag(Rx_signal_vp2)+tao/2,tao)-tao/2;
       Symbol_hat_vp = re + 1i*im;
       msg_hat_vp = demodulate(demod_obj,Symbol_hat_vp(:)/Scale);
       N_ebits_vp = N_ebits_vp + sum(msg_hat_vp ~= msg_bit);
    end
    BER(ii,1) = N_ebits_zf/N_tbits;
    BER(ii,2) = N_ebits_vp/N_tbits;
end
[beta_zf beta_vp;BER]
figure;semilogy(SNRdBs,BER(:,1),'-*');hold on
semilogy(SNRdBs,BER(:,2),'-ro');hold off;grid on


function [ out,arr,arr2] = precoding_vp( H,in,matt )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

W = H'*inv(H*H');
[NT,N_frame] = size(in);
tao = 2*sqrt(2);

out = zeros(NT,N_frame);
arr = zeros(1,N_frame);
arr2 = zeros(NT,N_frame);
test1 = zeros(NT,1);

for m=1:N_frame
    in_m = in(:,m);
    for n=1:6561
        test1 = W*(in_m+tao*matt(:,n));
        temp(n,1) = test1'*test1;
    end
    [min_val,min_loc] = min(temp);
    eff = in_m + tao*matt(:,min_loc);
    beta = real(sqrt(NT/(eff'*W'*W*eff)));
    arr(1,m) = beta;
    out(:,m) = beta*W*eff;
    arr2(:,m) = matt(:,min_loc);
end

end



function [out,beta ] = precoding_zf( H,in)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[NT,N_frame] = size(in);

W = H'*inv(H*H');
beta = sqrt(NT/trace(W*W'));
out = beta*W*in;

end

