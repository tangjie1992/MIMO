function[X_Estimated]= LRAD_MMSE2(H_complex,y,noise_variance,M,delta)
% Lattice Reduction with MMSE
% Input parameters
%     H_complex : complex channel matrix, nRxnT
%     y : complex received signal, nRx1
%     noise_variance : noise variance
% Output parameters
%     X_Estimated : estimated signal, nTx1

%MIMO-OFDM Wireless Communications with MATLAB㈢   Yong Soo Cho, Jaekwon Kim, Won Young Yang and Chung G. Kang
%?2010 John Wiley & Sons (Asia) Pte Ltd

alpha = sqrt(6/(M-1));
beta = sqrt(1.5/(M-1))*(sqrt(M)-1);
sigma = sqrt(noise_variance);

[Nr,Nt] = size(H_complex); N = 2*Nt;
N_frame = size(y,2);
Hr = real(H_complex); Hi = imag(H_complex);
H2 = [Hr -Hi; Hi Hr]; H3 = [H2;sigma*eye(N)]; % complex channel -> real channel  
y2 = [real(y);imag(y)]; y3 = [y2;zeros(N,N_frame)]+beta*H3*ones(2*Nr,N_frame); 
[Q,R,P,pp] = SQRD(H3);        % sorted QR decomposition 
[W,L,T] = original_LLL(Q,R,N,delta);  % W*L = Q*R*T

H3t = alpha*H3*P*T;             % H*P = Q*R
x3t = inv(H3t'*H3t)*H3t'*y3;  % MMSE detection
x3t2 = P*T*round(x3t);         % slicing
x3 = alpha*x3t2-beta*ones(2*Nr,N_frame);
X_Estimated = x3(1:Nt,:)+1i*x3(Nt+1:2*Nt,:); % real x -> complex x


function [Q,R,T] = original_LLL(Q,R,m,delta)
% Input parameters
%     Q : orthogonal matrix,  nRxnT
%     R : R with a large condition number
%     m : column length of H
%     delta : scaling variable
% Output parameters
%     Q : orthogonal matrix,  nRxnT
%     R : R with a small condition number
%     T : unimodular matrix
P=eye(m);  T=P;  k=2;
while (k <= m)
    for j = k-1:-1:1
        mu = round(R(j,k)/R(j,j));
        if (mu ~= 0)
            R(1:j,k)=R(1:j,k)-mu*R(1:j,j); 
            T(:,k)=T(:,k)-mu*T(:,j); 
        end
    end
    if (delta * R(k-1,k-1)^2 > R(k,k)^2 + R(k-1,k)^2)  % column change
       R(:,[k-1 k])=R(:,[k k-1]);
       T(:,[k-1 k])=T(:,[k k-1]);
%calculate Givens rotation matrix such that element R(k,k-1) becomes zero
       alpha = R(k-1,k-1)/sqrt(R(k-1:k,k-1).'*R(k-1:k,k-1));
       beta = R(k,k-1)/sqrt(R(k-1:k,k-1).'*R(k-1:k,k-1));
       theta = [alpha  beta; -beta  alpha];
       R(k-1:k,k-1:m)=theta*R(k-1:k,k-1:m);
       Q(:,k-1:k)=Q(:,k-1:k)*theta.';
       k=max([k-1 2]);
    else
       k=k+1;
    end
end

function [Q,R,P,p] = SQRD(H)
% Sorted QR decomposition
% Input parameter
%     H : complex channel matrix, nRxnT
% Output parameters
%     Q : orthogonal matrix, nRxnT
%     P : permutation matrix
%     p : ordering information
Nt=size(H,2);  Nr=size(H,1)-Nt;  R=zeros(Nt);
Q=H;   p=1:Nt;
for i=1:Nt 
    normes(i)=Q(:,i)'*Q(:,i); 
end
for i=1:Nt
   [mini,k_i]=min(normes(i:Nt)); k_i=k_i+i-1;
   R(:,[i k_i])=R(:,[k_i i]);
   p(:,[i k_i])=p(:,[k_i i]);
   normes(:,[i k_i])=normes(:,[k_i i]);
   Q(1:Nr+i-1,[i k_i])=Q(1:Nr+i-1,[k_i i]);
   % Wubben's algorithm: does not lead to
   % a true QR decomposition of the extended MMSE channel matrix
   % Q(Nr+1:Nr+Nt,:) is not triangular but permuted triangular
   R(i,i)=sqrt(normes(i));
   Q(:,i)=Q(:,i)/R(i,i);
   for k=i+1:Nt
      R(i,k)=Q(:,i)'*Q(:,k);
      Q(:,k)=Q(:,k)-R(i,k)*Q(:,i);
      normes(k)=normes(k)-R(i,k)*R(i,k)';
   end
end
P=zeros(Nt); 
for i=1:Nt
    P(p(i),i)=1;
end
