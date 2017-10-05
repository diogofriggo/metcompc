#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double a(double x0, double y0, double x1, double y1, double k);

int main()
{
    //char path[200] = "/Users/diogofriggo/Google Drive/UFRGS 8o Semestre/METODOS/metcompc/Aula9_0310/Results/ManySprings.txt";    
    int t, i, j, n = 2;
	double dt = 0.1, tmax = 200;
	double m1 = 1., m2 = 1., d = 1., k = 1., b = 1.;
	double vh_x[2], vh_y[2];
	double v_x[2], v_y[2];
	double x[2], y[2];
	for(i = 1; i < n-1; i++)
	{
		vh_x[i] = 0.;
		vh_y[i] = 0.;
		v_x[i] = 0.;
		v_y[i] = 0.;
		x[i] = 0.;
		y[i] = 0.;
	}
	//particle 0 starts at (0,b)
	//particle 1 starts at (10,0)
	y[0] = b;
	x[1] = 10.;
	v_x[0] = 1.;
    
    FILE *file = fopen("data.txt", "w");
    for(t = 0; t < tmax; t++)
    {
		//Report
        fprintf(file, "%.2f %.2f %.2f %.2f\n", x[0], y[0], x[1], y[1]);
        
        //Velocity-verlet
		vh_x[0] = v_x[0] + a(x[0], y[0], x[1], y[1], k)*dt/2.;
		vh_y[0] = v_y[0] + a(x[0], y[0], x[1], y[1], k)*dt/2.;
		vh_x[1] = v_x[1] + a(x[1], y[1], x[0], y[0], k)*dt/2.;
		vh_y[1] = v_y[1] + a(x[1], y[1], x[0], y[0], k)*dt/2.;
		//for(i = 0; i < n; i++)
		//{
        //    vh_x[i] = v_x[i] + a(x, y, k)*dt/2.;
		//	vh_y[i] = v_y[i] + a(x, y, k)*dt/2.;
		//}
        for(i = 0; i < n; i++)
		{
            x[i] = x[i] + vh_x[i]*dt;
            y[i] = y[i] + vh_y[i]*dt;
		}
        //for(i = 0; i < n; i++)
		//{
        //	v_x[i] = vh_x[i] + a(y, x, k)*dt/2.;
		//	v_y[i] = vh_y[i] + a(x, y, k)*dt/2.;
		//}
		v_x[0] = vh_x[0] + a(x[0], y[0], x[1], y[1], k)*dt/2.;
		v_y[0] = vh_y[0] + a(x[0], y[0], x[1], y[1], k)*dt/2.;
		v_x[1] = vh_x[1] + a(x[1], y[1], x[0], y[0], k)*dt/2.;
		v_y[1] = vh_y[1] + a(x[1], y[1], x[0], y[0], k)*dt/2.;
    }
    fclose(file);
    return 0;
}
//a_x, a_y
double a(double x0, double y0, double x1, double y1, double k)
{
	double xt = x0-x1;
	double yt = y0-y1;
	double r = sqrt(pow(xt,2.)+pow(yt,2.));
	return -k*r/pow(r,3.);
}




