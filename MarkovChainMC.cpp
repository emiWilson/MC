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

double EnergySq(){


	return Esq;
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

double MagnetizationSq(){
	
	return MagSq;

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
 	double dT = 0.1;
 	
 	double EqSteps = 10;

 	int numT = 4000;

 	double T[numT];

	for (int i = 0; i < numT; i++) 
	{
    	T[i] = maxT;
    	maxT = maxT - dT;
	}

	cout << endl;
	cout << "Initial Field" << endl;
	printSpinsToScreen();

	for (int t = 0; t < numT; t++){
		double temp = T[t];
		for (int eq = 0; eq < EqSteps; eq++){
			//NxM montecarlo steps per time step
			for (int n = 0; n < N; n++){
				for (int m = 0; m < M; m++){
					spinArr[n][m] = metropolis(temp, n, m);

				}
			}
		}
		double Eng = Energy() / (N*M);
		double Mag = (Magnetization()) / (N*M); //make indep if spins are all up or all  (ok to do? maybe not.)
		double Esq = EnergySq() /(N*M * N*M);
		double MagSq = MagnetizationSq()/(N*M * N*M);

		double Cv = (Esq - Eng*Eng) / (dT*dT) ; 
		double X = (MagSq - Mag*Mag) / dT ;

		Save(temp, Eng, Mag, Cv, X);

	}

	cout << "\n Final Field" << endl;
	printSpinsToScreen();

	
	return 0;
}