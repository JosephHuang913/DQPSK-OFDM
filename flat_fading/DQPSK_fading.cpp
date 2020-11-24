// Author: Chao-wang Huang
// Date: Friday, January 31, 2003
// BER Performance of DQPSK
// Jakes Rayleigh flat fading channel + AWGN noise

#define pi 3.141592653589793
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>

void AWGN_noise(float, double, double *);
void JakesFading(double, float, double, int, double *);
int Error_count(int, int);

const int sym_num = 10000000;
const float vc = 30.0; /* speed of vehicle in km/hr */
const double C = 3.0e8;  /* light speed */
const double fc = 2.0e9; /* carrier frequency */
const double sym_rate = 10.0e6;

int main(void)
{
	time_t  t, start, end;
	int  z, i, count, d_in[2], d_out[2];
   float Eb_No, err_rate, power_s, doppler;
   double ts, noise[2], signal_t[2], power_n, fade[2], phase_t, phase_r;
   double noise_power, tx_signal[2], rx_signal[2], phase_t_i, phase_r_i;
	FILE *ber, *records;

   start = time(NULL);
   printf("BER Performance for DQPSK under Rayleigh Flat Fading Channel\n");
   printf("This program is running. Don't close, please!\n\n");

	srand((unsigned) time(&t));
   doppler = (vc*1000.0/3600.0)/(C/fc);

	for(z=0; z<=35; z+=5)
   {
   	ts = 15.0;
		Eb_No = (float)z;
		noise_power = 0.5/pow(10.0,(double)Eb_No/10.0);
		count = 0;
		power_s = 0.0;
		power_n = 0.0;
      printf("Eb_No = %f, ", Eb_No);
      phase_t_i = 0.0;
      phase_r_i = 0.0;

		for(i=0; i<sym_num; i++)
   	{
         d_in[0] = random(2);
         d_in[1] = random(2);
         signal_t[0] = (1/sqrt(2.0))*(2*d_in[0]-1);
         signal_t[1] = (1/sqrt(2.0))*(2*d_in[1]-1);

         phase_t = phase_t_i + atan2(signal_t[1],signal_t[0]);
         phase_t_i = phase_t;
         tx_signal[0] = cos(phase_t);
         tx_signal[1] = sin(phase_t);

         power_s += (pow(tx_signal[0],2) + pow(tx_signal[1],2));

         JakesFading(fc, vc*1000/3600.0, ts, 2, &fade[0]);
         ts += 1.0/sym_rate;
         AWGN_noise(0.0, noise_power, &noise[0]);
         power_n += (pow(noise[0],2) + pow(noise[1],2));

         rx_signal[0] = tx_signal[0]*fade[0] - tx_signal[1]*fade[1] + noise[0];
         rx_signal[1] = tx_signal[0]*fade[1] + tx_signal[1]*fade[0] + noise[1];

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

      ber=fopen("ber_dqpsk_fading.log", "a");
	   records = fopen("records_dqpsk.log", "a");

      printf("Error rate = %e\n", err_rate);
      fprintf(ber, "%f %e\n", Eb_No, err_rate);
      fprintf(records, "DQPSK under flat fading channel with average Eb/N0 = %f dB\n", Eb_No);
      fprintf(records, "Speed of the vehicle = %f\n", vc);
	   fprintf(records, "Carrier Frequency = %e\n", fc);
      fprintf(records, "Maximum Doppler Frequency = %e\n", doppler);
      fprintf(records, "Transmission Symbol Rate = %e\n", sym_rate);
      fprintf(records, "Average signal power = %f\n", power_s);
	   fprintf(records, "Average noise power = %f\n", power_n);
   	fprintf(records, "Average bit error rate = %e\n", err_rate);
   	fprintf(records, "Symbol numbers of simulation = %d\n\n", sym_num);
      fclose(ber);
   	fclose(records);
	}

   end = time(NULL);
   records = fopen("records_dqpsk.log", "a");
   fprintf(records, "Total elapsed time: %.0f(sec)\n", difftime(end,start));
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

void JakesFading(double f_c/*Hz*/, float v/*m/s*/, double t/*s*/, int type, double *fade)
{
	const double C = 3.0e8;     // (m/s)
   const float Pi = 3.14159265358979;
   int n, N, N_o = 32;
   double lamda, w_m, beta_n, w_n, alpha, T_c2, T_s2, theta_n;

   lamda = C/f_c;     // wave length (meter)
   w_m = 2.0*Pi*v/lamda;    // maximum Doppler frequency
   N = 2*(2*N_o+1);

   switch(type)
   {
   	case 1:
   		alpha = 0.0;
         T_c2 = (double)N_o;
         T_s2 = (double)N_o + 1.0;
         break;
      case 2:
      	alpha = 0.0;
         T_c2 = (double)N_o + 1.0;
         T_s2 = (double)N_o;
         break;
      case 3:
      	alpha = Pi/4.0;
         T_c2 = (double)N_o + 0.5;
         T_s2 = (double)N_o + 0.5;
         break;
      default:
      	printf("\nInvalid type selection for Jake's fading channel model.\n");
         break;
   }

   if(v == 0.0)
   {
   	*(fade+0) = 1.0;
      *(fade+1) = 0.0;
   }
   else
   {
   	*(fade+0) = sqrt(1.0/T_c2)*cos(alpha)*cos(w_m*t);
      *(fade+1) = sqrt(1.0/T_s2)*sin(alpha)*cos(w_m*t);

      for(n = 1; n <= N_o; n++)
      {
      	switch(type)
         {
         	case 1:
            	beta_n = (double)n*Pi/((double)N_o+1.0);
               break;
            case 2:
            	beta_n = (double)n*Pi/(double)N_o;
               break;
            case 3:
            	beta_n = (double)n*Pi/(double)N_o;
               break;
         	default:
            	break;
         }
         w_n = w_m*cos(2.0*Pi*(double)n/(double)N);
//            theta_n = 2.0*Pi*((double)rand()/(double)RAND_MAX);  // random phase
			theta_n = 0.0;
         *(fade+0) += sqrt(2.0/T_c2)*cos(beta_n)*cos(w_n*t+theta_n);
         *(fade+1) += sqrt(2.0/T_s2)*sin(beta_n)*cos(w_n*t+theta_n);
		}
	}
}

int Error_count(int x, int y)
{
	if(x == y)
   	return 0;
   else
   	return 1;
}

