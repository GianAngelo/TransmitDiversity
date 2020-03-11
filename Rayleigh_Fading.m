%Gian Angelo Tria
%ECE335 
%PS4 Question 1 
clc;clear all; close all; 
%%
n = 12;
numbits = 2^n; %Number of bits being sent
data_0 = randi([0 1],numbits,1); %Random binary bits
data_1 = randi([0 1],numbits,1); 
fm = numbits/2;
SNR = linspace(0,50,51);

%BPSK Modulating the data using the functions
s_0 = transpose(BPSK(data_0));
s_1 = transpose(BPSK(data_1));


%Preallocating Arrays for number of errors and number
number_1 = zeros(5,length(SNR));
number_2 = zeros(5,length(SNR));
ratio_1 = zeros(5,length(SNR));
ratio_2 = zeros(5,length(SNR));
h = zeros(4,numbits);

%Making Rayleigh Fading
for i = 1:4
%Making gaussian distribution but negative frequencies are conjugated
%Shifted by fm to work with matlab matrix notation
gauss_0 = randn(1,2*fm) + 1i.*randn(1,2*fm);
gauss_0(1:fm) = conj(gauss_0(fm+1:end));
gauss_1 = randn(1,2*fm) + 1i.*randn(1,2*fm);
gauss_1(1:fm) = conj(gauss_1(fm+1:end));
%Computing doppler spectrum at baseband so fc = 0 
y = linspace(-fm,fm,numbits); 
dopp_spectrum = 1.5./ (pi*fm*sqrt(1-(y./fm).^2));
%Multiply them together
gauss_dopp_0 = gauss_0.*dopp_spectrum;
gauss_dopp_1 = gauss_1.*dopp_spectrum;
%The ends have a value of infinity so I put them to 0 instead
gauss_dopp_0([1 end]) = 0;
gauss_dopp_1([1 end]) = 0;

h(i,:) = sqrt(gauss_dopp_0.^2 + gauss_dopp_1.^2);
end

% 4 Different Rayleigh Fading Channels from the loop
h_0 = h(1,:);
h_1 = h(2,:);
h_2 = h(3,:);
h_3 = h(4,:);

%Preallocating for output of Maximum likelihood detector for each SNR
s_0_hat = zeros(numbits,length(SNR));
s_1_hat = zeros(numbits,length(SNR));

for j = 1:length(SNR)
    
%New Scheme 2Tx, 2Rx
%Recieved Signals through the different channels 
r_0 = awgn(h_0.*s_0 + h_1.*s_1,SNR(j),'measured');
r_1 = awgn(-h_0.*conj(s_1) + h_1.*conj(s_0),SNR(j),'measured');
r_2 = awgn(h_2.*s_0 + h_3.*s_1,SNR(j),'measured');
r_3 = awgn(-h_2.*conj(s_1) + h_3.*conj(s_0),SNR(j),'measured');

%Output of Combiner
s_tilda_0 = conj(h_0).*r_0 + h_1.*conj(r_1) + conj(h_2).*r_2 + h_3.*conj(r_3);

%Calulating distance for error
mag_err_stilda0_s0 = abs(s_tilda_0 - s_0);
mag_err_stilda0_s1 = abs(s_tilda_0 - s_1);

%Maximum Likelihood detector
for i = 1:numbits
    if mag_err_stilda0_s0(i) < mag_err_stilda0_s1(i)
        s_0_hat(i,j) = s_0(i);
    else
        s_0_hat(i,j) = s_1(i);
    end
end
%Demodulate and Compute BER
s_0_hat(:,j) = BPSKdemod(s_0_hat(:,j));
[number_1(5,j),ratio_1(5,j)] = biterr(s_0_hat(:,j),data_0);

%The steps for the every other scheme are pretty much the same, but if a
%channel is not present, that H is 0 so that term disappears. Also only
%computing errors for 1 signal to compare to the others with only
%transmitted signal.


%New Scheme 2Tx, 1Rx
r_0 = awgn(h_0.*s_0 + h_1.*s_1,SNR(j),'measured');
r_1 = awgn(-h_0.*conj(s_1) + h_1.*conj(s_0),SNR(j),'measured');

s_tilda_0 = conj(h_0).*r_0 + h_1.*conj(r_1);

mag_err_stilda0_s0 = abs(s_tilda_0 - s_0);
mag_err_stilda0_s1 = abs(s_tilda_0 - s_1);

for i = 1:numbits
    if mag_err_stilda0_s0(i) < mag_err_stilda0_s1(i)
        s_0_hat(i,j) = s_0(i);
    else
        s_0_hat(i,j) = s_1(i);
    end
end

s_0_hat(:,j) = BPSKdemod(s_0_hat(:,j));
[number_1(4,j),ratio_1(4,j)] = biterr(s_0_hat(:,j),data_0);


%Classical MRRC 1Tx 2Rx 
r_0 = awgn(h_0.*s_0,SNR(j),'measured');
r_1 = awgn(h_1.*s_0,SNR(j),'measured');

s_tilda_0 = conj(h_0).*r_0 + conj(h_1).*r_1;
mag_err_stilda0_0 = abs(s_tilda_0 - complex(1));
mag_err_stilda0_1 = abs(s_tilda_0 - complex(-1));

for i = 1:numbits
    if mag_err_stilda0_0(i) < mag_err_stilda0_1(i)
        s_0_hat(i,j) = complex(1);
    else
        s_0_hat(i,j) = complex(-1);
    end
end

s_0_hat(:,j) = BPSKdemod(s_0_hat(:,j));
[number_1(2,j),ratio_1(2,j)] = biterr(s_0_hat(:,j),data_0);


%Classical MRRC 1Tx 4Rx 
r_0 = awgn(h_0.*s_0,SNR(j),'measured');
r_1 = awgn(h_1.*s_0,SNR(j),'measured');
r_2 = awgn(h_2.*s_0,SNR(j),'measured');
r_3 = awgn(h_3.*s_0,SNR(j),'measured');

s_tilda_0 = conj(h_0).*r_0 + conj(h_1).*r_1 + conj(h_2).*r_2 + conj(h_3).*r_3;
mag_err_stilda0_0 = abs(s_tilda_0 - complex(1));
mag_err_stilda0_1= abs(s_tilda_0 - complex(-1));

for i = 1:numbits
    if mag_err_stilda0_0(i) < mag_err_stilda0_1(i)
        s_0_hat(i,j) = complex(1);
    else
        s_0_hat(i,j) = complex(-1);
    end
end
s_0_hat(:,j) = BPSKdemod(s_0_hat(:,j));
[number_1(3,j),ratio_1(3,j)] = biterr(s_0_hat(:,j),data_0);


% 1Tx 1Rx 

s_tilda_0 = conj(h_0).*r_0;
mag_err_stilda0_0 = abs(s_tilda_0 - complex(1));
mag_err_stilda0_1 = abs(s_tilda_0 - complex(-1));

for i = 1:numbits
    if mag_err_stilda0_0(i) < mag_err_stilda0_1(i)
        s_0_hat(i,j) = complex(1);
    else
        s_0_hat(i,j) = complex(-1);
    end
end
s_0_hat(:,j) = BPSKdemod(s_0_hat(:,j));
[number_1(1,j),ratio_1(1,j)] = biterr(s_0_hat(:,j),data_0);
end


semilogy(ratio_1(1,:))
hold on
semilogy(ratio_1(2,:))
hold on
semilogy(ratio_1(3,:))
hold on
semilogy(ratio_1(4,:))
hold on
semilogy(ratio_1(5,:))
axis([0 50 10^-6 1])
title('MRRC vs Two Branch Transmit Diversity through Rayleigh Fading')
legend('no diversity (1 Tx, 1Rx)','MRRC (1 Tx, 2Rx)','MRRC (1 Tx, 4Rx)','new scheme (2 Tx, 1 Rx)','new scheme (2 Tx, 2 Rx)');
xlabel('SNR')
ylabel('Pb bit error rate (BER)')

figure
scatterplot(s_0);
title('BPSK Modulated Signal')
xlabel('In Phase')
ylabel('Quadrature')

figure
scatterplot(h_0);
title('Example Rayleigh Channel')
xlabel('In Phase')
ylabel('Quadrature')

figure
plot(dopp_spectrum)
title('Doppler Spectrum for Rayleigh Fading')
xlabel('In Phase')
ylabel('Quadrature')

figure
scatterplot(gauss_0)
title('Example Gaussian with Negative Frequencies Conjugated')
xlabel('In Phase')
ylabel('Quadrature')
