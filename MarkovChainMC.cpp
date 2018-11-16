#include <iostream>
#include <fstream>
#include <math.h> /*pow*/
#include <random>

using namespace std;
FILE * pFile;
std::ofstream Data;
std::ifstream input;

int N; int M;
double maxT; double minT; int numT; double EqSteps; double dT;
double en; double mag; double en2; double mag2;
double Ene; double Mag;

int * s_ptr; double * t_ptr;

//Wipes the output for the two output files phaseField.txt and thermalField.txt
void restart(){

	Data.open("Data.txt", std::ofstream::out | std::ofstream::trunc);
	Data.close();

	srand(time(NULL)); //reset the random number generator
}
//read data in from a file
void getVariablesFromFile(){
  	cout <<"\nInput file: " << endl;
    string temporary;  // for para title
    ifstream inpf("Variables.txt");  // parameter file = para.in
    if (inpf.is_open()) {
    	//read N - first line of txt file
       	inpf >> temporary; inpf >> N;
        inpf >> temporary; inpf >> M;
        inpf >> temporary; inpf >> maxT;
        inpf >> temporary; inpf >> minT;
        inpf >> temporary; inpf >> numT;
	    inpf >> temporary; inpf >> EqSteps;
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

//Generate 1 or -1 randomly
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

//initialie the spin configuration
void initializeField(){
	for (int n = 0; n < N; n++){
		for (int m = 0; m < M; m++){
			* (s_ptr + (n * N) + m) = Rand();
		}
	}
}

//generate temperature array
void T_arr(double startT, double finalT, double steps){
	dT = (finalT - startT) / steps;
	for (int i = 0; i < steps; i++)
	{
    		* (t_ptr + i) = finalT;
    		finalT = finalT - dT; //decrement T by dT
	}
}

//the metropolis algorithm
int metropolis(double t, int n, int m){
	int s = spin(n, m);
	int NN = spin(BCx(n + 1), m) + spin(BCx(n - 1), m) + spin(n, BCy(m + 1)) + spin(n, BCy(m -1));
	int E_cost = 2*s*NN;

	double r = ((double) rand() / (RAND_MAX)); //random number between 0 and 1
	double beta = 1/t;

	//if flipping spin decreses energy then flip the spin
	if (E_cost < 0){
		s = - s;
	//else if special condtition is met flip spin
	}else if (r < exp(- E_cost * beta) ){
		s = -s;
	}
	return s;
}

//calculate Energy 
double Energy(){
	double E = 0;
	for (int n = 0; n < N; n ++){
		for (int m = 0; m < M; m++){
			int s = spin(n,m);
			//sum the spin of all nearest neighbours
			double NN = spin(BCx(n + 1), m) + spin(BCx(n - 1), m) + spin(n, BCy(m + 1)) + spin(n, BCy(m -1));
			//append energy of site to total energy E
			E = E - NN * s;
		}
	}
	return E;
}

//caculate magnetization by summming all the spins
double Magnetization(){
	double Mag = 0;
	for (int i = 0; i < N; i ++){
		for (int j = 0; j < M; j++){
			Mag = Mag + spin(i,j);
		}
	}
	return Mag;

}

//save data to .txt file
void Save(double temp, double E, double M, double Cv, double X){
	Data.open ("Data.txt", std::fstream::app);
	Data << temp << " " << E << " " << M << " " << Cv << " " << X << endl;
	Data.close();

}

//print spins config to console screen
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
	getVariablesFromFile(); //read in variables

	//initialize spin array
	int spinArr[N][M];
	s_ptr = & spinArr[0][0];
 	initializeField();

 	//print inital spin array to consol
	cout << endl;
	cout << "Initial Field" << endl;
	printSpinsToScreen();

	//relax system a bit before simulating
	for (int eq = 0; eq < EqSteps * 10; eq++){
			//NxM montecarlo steps per time step
			for (int n = 0; n < N; n++){
				for (int m = 0; m < M; m++){
					
					spinArr[n][m] = metropolis(maxT, n, m);
				}
			}
	}
	
	//Initialize T array  
	double T[numT];
	t_ptr = & T[0];
	T_arr(minT, maxT, numT);
 	
	//helpers for caluclating magnetic susceptability and heat capacity
	double n1 = 1.0 / (N * M * EqSteps);
	double n2 = 1.0 / (N * M * EqSteps * EqSteps);

	// equillibrate the system
	for (int t = 0; t < numT; t++){
		double temp = T[t];
		en = mag = en2 = mag2 = 0;

		double invT = 1.0 /temp; 
		double invT2 = invT * invT;


		//do EqSteps equillibration steps at each temperature step
		for (int eq = 0; eq < EqSteps; eq++){
			//NxM montecarlo steps per time step
			for (int n = 0; n < N; n++){
				for (int m = 0; m < M; m++){

					spinArr[n][m] = metropolis(temp, n, m);
				}
			}

			Ene = Energy(); //calc energy for equil step
			Mag = Magnetization(); //calc mag for equil step
			en = en + Ene; //append energy for equil step to energy for temp step
			mag = mag + Mag; //append mag for equil step to mag for time step
			
			en2 = en2 + (Ene * Ene); //calc energy squared for time step
			mag2 = mag2 + (Mag * Mag); //calc energy squared for time step
		}

		//Calculate heat capacity and magnetic susceptability via fluctuation dissipation thrm
		double Cv = (en2 * n1 - en * en * n2 ) * invT ; 
		double X = (mag2 * n1 - mag * mag * n2 ) * invT2;
		
		//save data for temp step to file
		Save(temp, Ene , Mag  , Cv, X);

	}

	//print final spin configuration to screen shot
	cout << "\n Final Field" << endl;
	printSpinsToScreen();

	return 0;
}