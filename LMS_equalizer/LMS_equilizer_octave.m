clc;
close all;
clear all;
TimeSlot=2e-3; %Transmit time duration
SNR = 18; %Signal to noise ratio
Rs = 185e3; % symbol rate
a=[1+0i 1/sqrt(2)+1i*1/sqrt(2) 1i -1/sqrt(2)+1i*1/sqrt(2) -1 -1/sqrt(2)-1i*1/sqrt(2) -1i 1/sqrt(2)-1i*1/sqrt(2) ];
Ak = a(randi(8,20000,1)); % 8 PSK sequence of 1000 samples
figure(1); plot(real(a),imag(a),'*'); title('8PSK constellation'); grid;

%% Channel creation and channel modelling
Rsym = Rs; M = 8;                  % Input symbol rate
Rbit = Rsym * log2(M);      % Input bit rate
Nos = 1;                    % Oversampling factor
ts = (1/Rsym) / Nos; 
pg=[0.9417 - 0.2214i -0.1596 - 0.0874i -0.0644 - 0.0163i -0.0645 - 0.0387i -0.0751 + 0.0467i];
pd=[0 2.0000e-06 4.0000e-06 6.0000e-06 8.0000e-06]./ts;
s=[Ak,zeros(1,1000)];
% N1=N2=800
%hence g is taken from 200 to 1800
for i=1:801
    y(i)=0;
    for n=200:1000
    g(n)=0;
    for k=1:5
        g(n)=g(n)+pg(k)*sinc(pd(k)-n+1000);
    end
        y(i)=y(i)+s(i-n+1000)*g(n);
    end
end
for j=802:20000
    y(j)=0;
    for n1=200:1800
    g(n1)=0;
    for k1=1:5
        g(n1)=g(n1)+pg(k1)*sinc(pd(k1)-n1+1000);
    end
        y(j)=y(j)+s(j-n1+1000)*g(n1);
    end
end
Rk=y
noise = (1/sqrt(2))*(randn(size(Rk)) + 1j*randn(size(Rk))); %Initial noise vector
P_s = var(Rk); % Signal power
P_n = var(noise); % Noise power
% Defining noise scaling factor based on the desired SNR:
noise_scaling_factor = sqrt(P_s/P_n./10.^(SNR./10)); 
Rk_noisy=Rk+noise*noise_scaling_factor; % Received signal
figure(2);  plot(real(Rk),imag(Rk),'*');legend('Received constellation')
hTap=11;%Channel Taps
beta = 0.001; % step-size of the algorithm
c_LMS = zeros(hTap,1); % equalizer coefficients, initializations
for i = (hTap+1)/2:length(Rk_noisy)-(hTap-1)/2 
rk = flipud(Rk_noisy(i-(hTap-1)/2 :i+(hTap-1)/2).'); % Received signal vector
Ek(i) = Ak(i) - c_LMS.'*rk; % Error signal, we assume a known symbol sequence
c_LMS = c_LMS + beta*Ek(i)*conj(rk); % LMS update !
end
%% MSE Equalizer%%%%%%%%%%%%%%%%%%%%%%%LMS Algorithm%%%%%%%%%%%%%%%%
% Initialization
FII = zeros(hTap,hTap); % Autocorrelation Matrix initialization
alfa = zeros(hTap,1); % Cross-correlation vector initialization
% Estimating FII and alfa using sample estimates based on training symbols
% Notice that here we use all the generated data as training symbols
for i = (hTap+1)/2:length(Rk_noisy)-(hTap-1)/2 %16:length(Rk_noisy)-15, 
% rk = flipud(Rk_noisy(i-15:i+15).'); % Received signal vector
rk = flipud(Rk_noisy(i-(hTap-1)/2 :i+(hTap-1)/2).'); % Received signal vector
FII = FII + rk*conj(rk).'; % Autocorrelation matrix
alfa = alfa + Ak(i)*conj(rk); 
end
% FII
FII = FII/(length(Rk_noisy)-(hTap-1)); % Final sample estimate of the autocorrelation matrix
alfa = alfa/(length(Rk_noisy)-(hTap-1)); % Final sample estimate of the cross-correlation vector
c_MSE = inv(conj(FII))*alfa; % Equalizer coefficients
%%
% Plotting
figure(3);
subplot(121);
y_LMS=filter(c_LMS,1,Rk);
plot(y_LMS,'rx'); 
title('LMS Equalized Constellation');
grid on;
subplot(122);
y_MSE=filter(c_MSE,1,Rk);
plot(y_MSE((hTap+1)/2:end),'rx');
title('MSE Equalized Constellation');
grid on;
