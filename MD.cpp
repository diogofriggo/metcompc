#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
using namespace std;

double rd(){ // random number between 0 and 1
	int a=rand()%100;
	double rd=0.01*a;
	
	return rd;
}

double MaxBol(int j , double T){ // Maxwell-Boltzmann distribution
	int x=0, y=1, z=2;
	double u1=rd();
	double u2=rd();
	double u3=rd();
	double v;
	
	if(j==x){
		v=sqrt(-2.0*T*log(1-u1))*cos(2.0*M_PI*u2)*sin(2.0*M_PI*u3);
	}
	if(j==y){
		v=sqrt(-2.0*T*log(1-u1))*sin(2.0*M_PI*u2)*sin(2.0*M_PI*u3);
	}
	if(j==z){
		v=sqrt(-2.0*T*log(1-u1))*cos(2.0*M_PI*u3);
	}
	
	return v;
}

double magnitude(double x , double y , double z){ // calculating magnitude
	double magnitude=sqrt(x*x + y*y + z*z);
	
	return magnitude;
}

int main(int argc, char* argv[]){
	double dt=0.005; // time-step
	double tE=500;   // end-time
	
	// System:
	int N=atoi(argv[1]);    // particles
	double T=atof(argv[2]); // temperature
	double L=25;            // square system
	int x=0 , y=1 , z=2;
	double r[N][3];         // location[particle][coordinate]
	double v[N][3];         // velocity[particle][coordinate]
	double vh[N][3];        // velocity[particle][coordinate] at 0.5*dt
	double a[N][3];         // acceleration[particle][coordinate]
	double m=1.0;           // mass
	
	// Force:
	double d;                   // distance between particles
	double dmin=0.85;           // minimum force range for stability
	double dmax=sqrt(2)*L;              // maximum force range for computation time
	double epsilon=1.0;
	double sigma=1/1.12246;
	double c1=-24.0*epsilon/m;
	double sigma6=pow(sigma,6);
	double d6 , d13i;
	double c2;
	double aik;
	
	
	// Initial Conditions:
	double dx=1.0*L/(1.0*sqrt(N));
	double dy=1.0*L/(1.0*sqrt(N));
	int h=0;
	for(int i=0 ; i<=sqrt(N)-1 ; i++){
		for(int j=0 ; j<=sqrt(N)-1 ; j++){
			r[h][x]=i*dx;
			r[h][y]=j*dy;
			r[h][z]=0.0;
			h++;
		}
	}
	
	for(int i=0 ; i<=N-1 ; i++){
		v[i][x]=MaxBol(x,T);
		v[i][y]=0.0;
		v[i][z]=0.0;
	}
	
	for(int i=0 ; i<=N-1 ; i++){
		for(int j=x ; j<=z ; j++){
			a[i][j]=0.0;
		}
	}	
	for(int i=0 ; i<=N-1 ; i++){
		for(int k=i+1 ; k<=N-1 ; k++){
			d=magnitude(r[i][x]-r[k][x] , r[i][y]-r[k][y] , r[i][z]-r[k][z]);
			if(d>=dmin && d<=dmax){
				d6=(d*d*d*d*d*d);
				d13i=1/(d6*d6*d);
				c2=c1*sigma6*d13i*(d6-2.0*sigma6);
				for(int j=x ; j<=z ; j++){
					aik=c2*(r[i][j]-r[k][j]);
					a[i][j]=a[i][j]+aik;
					a[k][j]=a[k][j]-aik;
				}
			}
		}
	}
		
	for(int i=0 ; i<=N-1 ; i++){
		cout << r[i][x] << "\t" << r[i][y] << "\t" << r[i][z] << "\t" << magnitude(v[i][x],v[i][y],v[i][z]) << "\t";
	}
	cout << "\n";

	// Time-evolution:
	int Nt=tE/dt; 
	int out=1/dt;
	int p=5;
	for(int t=1 ; t<=Nt ; t++){
		for(int i=0 ; i<=N-1 ; i++){
			for(int j=x ; j<=z ; j++){
				vh[i][j]=v[i][j] + 0.5*dt*a[i][j];
			}
			
			for(int j=x ; j<=z ; j++){
				r[i][j]=r[i][j]+vh[i][j]*dt;
				if(r[i][j]<0){r[i][j]=L+r[i][j];} // periodic BC
				if(r[i][j]>L){r[i][j]=r[i][j]-L;} // periodic BC
			}
		}	
		
		for(int i=0 ; i<=N-1 ; i++){
			for(int j=x ; j<=z ; j++){
				a[i][j]=0.0;
			}
		}	
		for(int i=0 ; i<=N-1 ; i++){
			for(int k=i+1 ; k<=N-1 ; k++){
				d=magnitude(r[i][x]-r[k][x] , r[i][y]-r[k][y] , r[i][z]-r[k][z]);
				if(d>=dmin && d<=dmax){
					d6=(d*d*d*d*d*d);
					d13i=1/(d6*d6*d);
					c2=c1*sigma6*d13i*(d6-2.0*sigma6);
					for(int j=x ; j<=z ; j++){
						aik=c2*(r[i][j]-r[k][j]);
						a[i][j]=a[i][j]+aik;
						a[k][j]=a[k][j]-aik;
					}
				}
			}
		}
		
		for(int i=0 ; i<=N-1 ; i++){	
			for(int j=x ; j<=z ; j++){
				v[i][j]=vh[i][j] + 0.5*dt*a[i][j];
			}
		}
		
		if(t==out){
			for(int i=0 ; i<=N-1 ; i++){
				cout << r[i][x] << "\t" << r[i][y] << "\t" << r[i][z] << "\t" << magnitude(v[i][x],v[i][y],v[i][z]) << "\t";
			}
			cout << "\n";
			
			out=out+1/0.005;
		}
		
		if((100*t/Nt)==p){
			cerr << 100*t/Nt << "% \n";
			p=p+5;
		}
	}	
}