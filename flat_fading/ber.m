close;
clear;
ber1 = fopen('ber_dqpsk_ofdm_fading_256.log', 'r');
ber = fscanf(ber1, '%e');
i=1:1:8;
Eb_No(i) = ber(2*i-1);
err_rate(i) = ber(2*i);
semilogy(Eb_No, err_rate, '-r^');
hold on;
grid on;
title('BER Performance of DQPSK + OFDM under Rayleigh Flat Fading Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
fclose(ber1);

ber1 = fopen('ber_dqpsk_ofdm_fading_128.log', 'r');
ber = fscanf(ber1, '%e');
i=1:1:8;
Eb_No(i) = ber(2*i-1);
err_rate(i) = ber(2*i);
semilogy(Eb_No, err_rate, '-bs');
fclose(ber1);

ber1 = fopen('ber_dqpsk_fading.log', 'r');
ber = fscanf(ber1, '%e');
i=1:1:8;
Eb_No(i) = ber(2*i-1);
err_rate(i) = ber(2*i);
semilogy(Eb_No, err_rate, '-k.');
fclose(ber1);

legend('DQPSK + OFDM, FFT: 256', 'DQPSK + OFDM, FFT: 128', 'DQPSK', 'f_c = 2.0 GHz', 'f_d = 55.56 Hz', 'Rs = 10 MHz');
%print -djpeg100 ber_dqpsk_ofdm_fading.jpg
