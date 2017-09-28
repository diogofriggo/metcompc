#include <stdio.h>
#include <math.h>

double aon1(double M, double k, double k12, double x1, double x2);
double aon2(double M, double k, double k12, double x1, double x2);

void main()
{
	double dt = 0.1;
	double M = 1.;

	double x1 = 0.3;
	double x1p = 0.0;
	double v1 = 0.1;
	double v1p = 0.0;
	
	double x2 = 0.6;
	double x2p = 0.0;
	double v2 = -0.1;
	double v2p = 0.0;
	
	double k = 1.;
	double k12 = 1.;
	int t = 0.0;
	FILE *file = fopen("data.txt", "w");
	fprintf(file, "%.2f %.2f %.2f %.2f %.2f\n", t*dt, x1, v1, x2, v2);
	for(t = 0; t < 1000; t++)
	{
		v1p = v1 + a(M, k, k12, x1, x2)*dt/2.;
		xp = x + v*dt;
		vp = v + a(M, k, k12, x1, x2)*dt/2.;

		v = v + a(M, k, k12, x1, x2)*dt/2.;
		xp = x + v*dt;
		vp = v + a(M, k, k12, x1, x2)*dt/2.;

		x1 = x1p;
		v1 = v1p;
		x2 = x2p;
		v2 = v2p;
		fprintf(file, "%.2f %.2f %.2f %.2f %.2f\n", t*dt, x1, v1, x2, v2);
	}
	fclose(file);
}

double aon1(double M, double k, double k12, double x1, double x2)
{
	return -k*x1-k12*(x1-x2)/M;
}

double aon2(double M, double k, double k12, double x1, double x2)
{
	return -k*x2-k12*(x2-x1)/M;
}

//		x1p = x1 + v1*dt + fon1(k, k12, x1, x2)*pow(dt,2.)/2.;
//		x2p = x2 + v2*dt + fon2(k, k12, x1, x2)*pow(dt,2.)/2.;
		
//		v1p = v1 + (fon1(k, k12, x1, x2)+fon1(k, k12, x1p, x2p))*dt/2.;
//		v2p = v2 + (fon2(k, k12, x1, x2)+fon2(k, k12, x1p, x2p))*dt/2.;

//		x1 = x1p;
//		v1 = v1p;
//		x2 = x2p;
//		v2 = v2p;
//		fprintf(file, "%.2f %.2f %.2f %.2f %.2f\n", t*dt, x1, v1, x2, v2);

