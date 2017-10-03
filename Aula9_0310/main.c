#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double a(double *x, double k_left, double k_right, int index, double M, int n);

int main()
{
    double dt = 0.01;
    double tmax = 10000;
    double M = 1.;
    double L = 1.;
    double k = 1.;
    double k01 = 1.;
    int t = 0;
	int i = 0;
    const int n = 3;
    
    double x[n];
    double v[n];
    double vh[n];
    
    x[0] = 0.1;
    v[0] = 0.;
    vh[0] = 0.;
    
    x[1] = 0.;
    v[1] = 0.;
    vh[1] = 0.;

	x[2] = 0.;
    v[2] = 0.;
    vh[2] = 0.;
    
    FILE *file = fopen("data.txt", "w");
    fprintf(file, "%.2f %.2f %.2f %.2f %.2f\n", t*dt, x[0]+d, v[0], x[1]+L/.2, v[1]);
    for(t += dt; t < tmax; t++)
    {
		for(i = 0; i < n; i++)
			vh[i] = v[i] + a(x, k_left, k_right, i, M, n)*dt/2.;
		for(i = 0; i < n; i++)
			x[i] = x[i] + vh[i]*dt;
		for(i = 0; i < n; i++)
			v[i] = vh[i] + a(x, k_left, k_right, i, M, n)*dt/2.;
		fprintf(file, "%.2f ", t*dt);
		for(i = 0; i < n; i++)
	        fprintf(file, "%.2f %.2f ", x[i]+i*d, v[i]);
		fprintf(file, "\n");
    }
    fclose(file);
    return 0;
}

double a(double *x, double k_left, double k_right, int middle, double M, int n, double d)
{
	int left = middle-1;
	int right = middle+1;
	if(index == 0) //x_left fixed
		return k_right*(x[right]+right*d) - (k_left+k_right)*(x[middle]+middle*d);
	if(index == n) //x_right fixed
		return k_left*(x[left]+left*d) - (k_left+k_right)*(x[middle]+middle*d);
	return k_left*(x[left]+left*d) + k_right*(x[right]+right*d) - (k_left+k_right)*(x[middle]+middle*d);
}
