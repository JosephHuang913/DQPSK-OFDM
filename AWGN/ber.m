close;
clear;
ber1 = fopen('ber_dqpsk_ofdm_awgn.log', 'r');
ber_awgn = fscanf(ber1, '%e');
i=1:1:16;
Eb_No(i) = ber_awgn(2*i-1);
err_rate(i) = ber_awgn(2*i);
semilogy(Eb_No, err_rate, '-r*');
hold on;
grid on;
%semilogy(Eb_No, err_rate, 'r*');
title('BER Performance of DQPSK + OFDM under AWGN Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
fclose(ber1);

ber1 = fopen('ber_dqpsk_awgn.log', 'r');
ber_awgn = fscanf(ber1, '%e');
i=1:1:16;
Eb_No(i) = ber_awgn(2*i-1);
err_rate(i) = ber_awgn(2*i);
semilogy(Eb_No, err_rate, '-k.');
fclose(ber1);

legend('DQPSK + OFDM', 'DQPSK', 'FFT: 128');
%print -djpeg100 ber_dqpsk_ofdm_awgn.jpg;
