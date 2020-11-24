// Author: Chao-wang Huang
// Date: Friday, January 31, 2003
// BER Performance of DQPSK + OFDM
// Jakes Rayleigh flat fading channel + AWGN noise

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define pi 3.141592653589793
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>

void dfour1(double *, unsigned long, int);
void AWGN_noise(float, double, double *);
void JakesFading(double, float, double, int, double *);
int Error_count(int, int);

const int n_fft = 256;
const int ofdm_sym_num = 40000;
const float vc = 30.0; /* speed of vehicle in km/hr */
const double C = 3.0e8;  /* light speed */
const double fc = 2.0e9; /* carrier frequency */
const double sym_rate = 10.0e6;

int main(void)
{
	time_t  t, start, end;
	int  z, i, j, total_samples, count, d_in[512], d_out[512];
   float Eb_No, err_rate, power_s, doppler;
   double ts, noise[2], signal_t[512], power_n, fade[2], phase_t[256], phase_r[256];
   double noise_power, tx_signal[513], rx_signal[513], phase_t_i[256], phase_r_i[256];
	FILE *ber, *records;

   start = time(NULL);
   printf("BER Performance for DQPSK + OFDM under Rayleigh Flat Fading Channel\n");
   printf("This program is running. Don't close, please!\n\n");

	srand((unsigned) time(&t));
   total_samples = 2*n_fft*ofdm_sym_num;
   doppler = (vc*1000.0/3600.0)/(C/fc);

	for(z=0; z<=35; z+=5)
   {
   	ts = 15.0;
		Eb_No = (float)z;
		noise_power = 0.5/(float)n_fft/pow(10.0,(double)Eb_No/10.0);
		count = 0;
		power_s = 0.0;
		power_n = 0.0;
      printf("Eb_No = %f, ", Eb_No);

      for(j=0; j<n_fft; j++)
      {
      	phase_t_i[j] = 0.0;
         phase_r_i[j] = 0.0;
      }

		for(i=0; i<ofdm_sym_num; i++)
   	{
      	for(j=0; j<n_fft; j++)
         {
         	d_in[2*j] = random(2);
            d_in[2*j+1] = random(2);
            signal_t[2*j] = (1/sqrt(2.0))*(2*d_in[2*j]-1);
            signal_t[2*j+1] = (1/sqrt(2.0))*(2*d_in[2*j+1]-1);

            phase_t[j] = phase_t_i[j] + atan2(signal_t[2*j+1],signal_t[2*j]);
            phase_t_i[j] = phase_t[j];

            tx_signal[2*j+1] = cos(phase_t[j]);
            tx_signal[2*j+2] = sin(phase_t[j]);
         }

         dfour1(&tx_signal[0],256,1);  // IFFT for transmitted signal must be multiplied by n_fft^(-1)

         for(j=0; j<n_fft; j++)
         {
         	tx_signal[2*j+1] = (tx_signal[2*j+1])/(float)n_fft;
            tx_signal[2*j+2] = (tx_signal[2*j+2])/(float)n_fft;
            power_s += (pow(tx_signal[2*j+1],2) + pow(tx_signal[2*j+2],2));

            JakesFading(fc, vc*1000/3600.0, ts, 2, &fade[0]);
            ts += 1.0/sym_rate;
            AWGN_noise(0.0, noise_power, &noise[0]);
            power_n += (pow(noise[0],2) + pow(noise[1],2));

            rx_signal[2*j+1] = tx_signal[2*j+1]*fade[0] - tx_signal[2*j+2]*fade[1] + noise[0];
            rx_signal[2*j+2] = tx_signal[2*j+1]*fade[1] + tx_signal[2*j+2]*fade[0] + noise[1];
         }

         dfour1(&rx_signal[0],256,-1);

         for(j=0; j<n_fft; j++)
         {
         	phase_r[j] = atan2(rx_signal[2*j+2],rx_signal[2*j+1]) - phase_r_i[j];
            phase_r_i[j] = atan2(rx_signal[2*j+2],rx_signal[2*j+1]);

            if(cos(phase_r[j])>=0.0)
            	d_out[2*j] = 1;
            else
            	d_out[2*j] = 0;

            if(sin(phase_r[j])>=0.0)
            	d_out[2*j+1] = 1;
            else
            	d_out[2*j+1] = 0;

            count += Error_count(d_in[2*j], d_out[2*j]);
            count += Error_count(d_in[2*j+1], d_out[2*j+1]);
         }
		}
   	err_rate = count/(float)total_samples;
      power_s = 2*power_s/(float)total_samples;
   	power_n = 2*power_n/(double)total_samples;

      ber=fopen("ber_dqpsk_ofdm_fading_256.log","a");
   	records = fopen("records_dqpsk_ofdm_fading_256.log", "a");

      printf("Error rate = %e\n", err_rate);
      fprintf(ber, "%f %e\n", Eb_No, err_rate);
      fprintf(records, "DQPSK + OFDM under flat fading channel\n");
      fprintf(records, "Average Eb/N0 = %f dB\n", Eb_No);
      fprintf(records, "Speed of the vehicle = %f\n", vc);
	   fprintf(records, "Carrier Frequency = %e\n", fc);
      fprintf(records, "Maximum Doppler Frequency = %e\n", doppler);
      fprintf(records, "Transmission Symbol Rate = %e\n", sym_rate);
      fprintf(records, "Average signal power = %f\n", power_s);
	   fprintf(records, "Average noise power = %f\n", power_n);
   	fprintf(records, "Average bit error rate = %e\n", err_rate);
   	fprintf(records, "Sample numbers of simulation = %d\n", total_samples);
      end = time(NULL);
      fprintf(records, "Total time elapsed: %.0f(sec)\n\n", difftime(end, start));
		fclose(ber);
   	fclose(records);
	}

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

void dfour1(double data[],unsigned long nn,int isign)
//double data[];
//unsigned long nn;
//int isign;
{
	unsigned long n,mmax,m,j,istep,i;
   double wtemp,wr,wpr,wpi,wi,theta;
   double tempr,tempi;

   n=nn << 1;
   j=1;
   for (i=1;i<n;i+=2)
   {
   	if (j > i)
      {
      	SWAP(data[j],data[i]);
         SWAP(data[j+1],data[i+1]);
      }
      m=n >> 1;
      while (m >= 2 && j > m)
      {
      	j -= m;
         m >>= 1;
      }
      j += m;
   }
   mmax=2;
   while (n > mmax)
   {
   	istep=mmax << 1;
      theta=isign*(6.28318530717959/mmax);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (m=1;m<mmax;m+=2)
      {
      	for (i=m;i<=n;i+=istep)
         {
         	j=i+mmax;
            tempr=wr*data[j]-wi*data[j+1];
            tempi=wr*data[j+1]+wi*data[j];
            data[j]=data[i]-tempr;
            data[j+1]=data[i+1]-tempi;
            data[i] += tempr;
            data[i+1] += tempi;
         }
         wr=(wtemp=wr)*wpr-wi*wpi+wr;
         wi=wi*wpr+wtemp*wpi+wi;
      }
      mmax=istep;
   }
}
#undef SWAP

