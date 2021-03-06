// Author: Chao-wang Huang
// Date: Saturday, February 1, 2003
// BER Performance of DQPSK under AWGN channel

#define pi 3.141592653589793
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>

void AWGN_noise(float, double, double *);
int Error_count(int, int);

const int sym_num = 1000000;

int main(void)
{
	time_t  t, start, end;
	int  z, i, count, d_in[2], d_out[2];
   float Eb_No, err_rate, power_s;
   double noise[2], signal_t[2], power_n, phase_t, phase_r;
   double noise_power, tx_signal[2], rx_signal[2], phase_t_i, phase_r_i;
	FILE *ber, *records;

   start = time(NULL);
   printf("BER Performance for DQPSK under AWGN Channel\n");
   printf("This program is running. Don't close, please!\n\n");

	srand((unsigned) time(&t));

	for(z=0; z<=15; z++)
   {
		Eb_No = (float)z;
		noise_power = 0.5/pow(10.0,(double)Eb_No/10.0);
		count = 0;
		power_s = 0.0;
		power_n = 0.0;
      printf("Eb_No = %f, ", Eb_No);
      phase_t = 0.0;
      phase_r_i = 0.0;

		for(i=0; i<sym_num; i++)
   	{
         d_in[0] = random(2);
         d_in[1] = random(2);
         signal_t[0] = (1/sqrt(2.0))*(2*d_in[0]-1);
         signal_t[1] = (1/sqrt(2.0))*(2*d_in[1]-1);

         phase_t_i = phase_t;
         phase_t = phase_t_i + atan2(signal_t[1],signal_t[0]);
         tx_signal[0] = cos(phase_t);
         tx_signal[1] = sin(phase_t);

         power_s += (pow(tx_signal[0],2) + pow(tx_signal[1],2));

         AWGN_noise(0.0, noise_power, &noise[0]);
         power_n += (pow(noise[0],2) + pow(noise[1],2));

         rx_signal[0] = tx_signal[0] + noise[0];
         rx_signal[1] = tx_signal[1] + noise[1];

         phase_r = atan2(rx_signal[1],rx_signal[0]) - phase_r_i;
         phase_r_i = atan2(rx_signal[1],rx_signal[0]);

         if(cos(phase_r)>=0.0)
         	d_out[0] = 1;
         else
         	d_out[0] = 0;

         if(sin(phase_r)>=0.0)
         	d_out[1] = 1;
         else
         	d_out[1] = 0;

         count += Error_count(d_in[0], d_out[0]);
         count += Error_count(d_in[1], d_out[1]);
		}
   	err_rate = count/(float)(2*sym_num);
      power_s = power_s/(float)sym_num;
   	power_n = power_n/(double)sym_num;

      ber=fopen("ber_dqpsk_awgn.log","a");
   	records = fopen("records_dqpsk.log", "a");

      printf("Error rate = %e\n", err_rate);
      fprintf(ber, "%f %e\n", Eb_No, err_rate);
      fprintf(records, "DQPSK under AWGN channel\n");
      fprintf(records, "Eb/N0 = %f dB\n", Eb_No);
      fprintf(records, "Average signal power = %f\n", power_s);
	   fprintf(records, "Average noise power = %f\n", power_n);
   	fprintf(records, "Average bit error rate = %e\n", err_rate);
   	fprintf(records, "Symbol numbers of simulation = %d\n\n", sym_num);
		fclose(ber);
   	fclose(records);
	}

   records = fopen("records_dqpsk.log", "a");
   end = time(NULL);
   fprintf(records, "Totle elapsed time: %.0f(sec)\n", difftime(end,start));
   fclose(records);
	printf("This program is ended. Press any key to continue.\n");
   getch();

   return 0;
}

void AWGN_noise(float mu, double variance, double *noise)
{
	const float Pi = 3.14159265358979;
   double u1, u2;
   do
   {
   	u1 = (double)rand()/(double)RAND_MAX;
      u2 = (double)rand()/(double)RAND_MAX;
   }
   while(u1 == 0.0 || u2 == 0.0);

   *(noise+0) = (sqrt(-2.0*log(u1))*cos(2*Pi*u2))*sqrt(variance/2.0)+mu/sqrt(2.0);
   *(noise+1) = (sqrt(-2.0*log(u1))*sin(2*Pi*u2))*sqrt(variance/2.0)+mu/sqrt(2.0);
}

int Error_count(int x, int y)
{
	if(x == y)
   	return 0;
   else
   	return 1;
}

