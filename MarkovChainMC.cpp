#include <iostream>
#include <fstream>
#include <math.h> /*pow*/
#include <random>

using namespace std;
FILE * pFile;

std::ofstream Data;

std::ifstream input; //if wanted to read file in.

int N = 30;
int M = N; //square lattice

int * s_ptr;

//Wipes the output for the two output files phaseField.txt and thermalField.txt
void restart(){

	Data.open("Data.txt", std::ofstream::out | std::ofstream::trunc);
	Data.close();

	srand(time(NULL)); //reset the random number generator
}

void getVariablesFromFile(){


  	cout <<"\nInput file: " << endl;

    string temp;  // for para title
    ifstream inpf("variables.txt");  // parameter file = para.in
    if (inpf.is_open()) {
    	//read N - first line of txt file
       	inpf >> temp; inpf >> N;
        M = N;
 

        inpf.close();

    } else {

        cout << "Failed to locate parameter file (\"para.in\")" << endl;

    }
    
    
}

/*** Periodic Boundary Conditons ***/
int BCx(int n){
	return (n + N) % N ;
}

int BCy(int m){
	return (m + M) % N ;
}

/*** Generate 1 or -1 randomly ***/
int Rand(){
	int r = rand() % 2;
	if (r == 0)
		r = -1;

	return r;
}

//gets values in the phi array
int spin(int i, int j){

	int element = * (s_ptr + (BCx(i) * N) + BCy(j));
	
	return element;
}

/*** Initialie the spin configuration ***/
void initializeField(){

	for (int n = 0; n < N; n++){
		for (int m = 0; m < M; m++){
			* (s_ptr + (n * N) + m) = Rand();
		}
	}

}

int metropolis(double t, int n, int m){
	int s = spin(n, m);
	int NN = spin(BCx(n + 1), m) + spin(BCx(n - 1), m) + spin(n, BCy(m + 1)) + spin(n, BCy(m -1));
	int E_cost = 2*s*NN;

	double r = ((double) rand() / (RAND_MAX)); //random number between 0 and 1
	double beta = 1/t;

	if (E_cost < 0){
		s = - s;

	}else if (r < exp(- E_cost * beta) ){
		s = -s;

	}
	return s;
}

double Energy(){
	double E = 0;
	for (int n = 0; n < N; n ++){
		for (int m = 0; m < M; m++){
			int s = spin(n,m);
			double NN = spin(BCx(n + 1), m) + spin(BCx(n - 1), m) + spin(n, BCy(m + 1)) + spin(n, BCy(m -1));
			E = E - NN*s;
		}
	}

	return E;
}

double Magnetization(){
	double Mag = 0;
	for (int i = 0; i < N; i ++){
		for (int j = 0; j < M; j++){
			Mag = Mag + spin(i,j);
		}
	}
	return Mag;

}

void Save(double temp, double E, double M, double Cv, double X){
	
	Data.open ("Data.txt", std::fstream::app);
	Data << temp << " " << E << " " << M << " " << Cv << " " << X << endl;
	Data << endl;
	Data.close();

}


void printSpinsToScreen(){
	for (int i = 0; i < N ; i ++){
		for (int j = 0; j < M; j ++){
			int s= spin(i,j);
			if (s == 1)
				cout << "+"<<s << " ";
			else
				cout << s << " ";

		}
		cout << endl;
	}

	cout << endl;
}
int main(){

	restart(); //make sure starting with empty output files

	int spinArr[N][M];
	s_ptr = & spinArr[0][0];
 	initializeField();


 	/*** Initialize T array ***/ 
 	double maxT = 4;

 	double EqSteps = 1000;

 	int numT = 100;
 	double dT = maxT / numT;

 	double T[numT];

	for (int i = 0; i < numT; i++) 
	{
    	T[i] = maxT;
    	maxT = maxT - dT;
	}

	cout << endl;
	cout << "Initial Field" << endl;
	printSpinsToScreen();

	double en;
	double mag;
	double en2;
	double mag2;

	double n1 = 1.0 / (N*M*EqSteps);
	double n2 = 1.0 / (N*M*EqSteps*EqSteps);

	for (int t = 0; t < numT; t++){
		double temp = T[t];
		en = mag = en2 = mag2 = 0;

		double invT = 1.0/temp; 
		double invT2 = invT * invT;

		//do EqSteps equillibration steps at each temperature step
		for (int eq = 0; eq < EqSteps; eq++){
			//NxM montecarlo steps per time step
			for (int n = 0; n < N; n++){
				for (int m = 0; m < M; m++){
					
					spinArr[n][m] = metropolis(temp, n, m);

				}
			}
			double Ene = Energy();
			double Mag = Magnetization();

			en = en + Ene;
			mag = mag + Mag;
			en2 = en2 + Ene*Ene;
			mag2 = mag2 + Mag * Mag;

		}


		double Eng = en * n1;
		double Mag = mag * n1;

		double Cv = (en2*n1 - en * en * n2) / (invT2) ; 
		double X = (mag2*n1 - mag * mag * n2) / (invT);
		//cout << temp << endl;
		Save(temp, Eng, Mag, Cv, X);

	}

	cout << "\n Final Field" << endl;
	printSpinsToScreen();

	
	return 0;
}