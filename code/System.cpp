#include "System.h"
#include "Type.h"
#include <iostream>
#include <fstream> // zapis do pliku
#include <cmath>

using namespace std;

System::System(double w, double h, int x, int y, double gx, double gy, double v):
	width(w), height(h), nx(x), ny(y), gx(gx), gy(gy), vis(v){

		dx = width/nx;
		dy = height/ny;

		// tablice

		U1.init(nx, ny);
		V1.init(nx, ny);
		P1.init(nx, ny);
		D1.init(nx, ny);
		Z1.init(nx, ny);
		F1.init(nx, ny);


		// dokladnosc w relaksacji cisnienia
		EpsP = 0.0002; 

	 	DiffPMax = 2.0;

		neighSum = 0;
		
		iters = 0;

		// do zmiany w razie potrzeby
		frameTime = 1.0/2.0; 

		//poczatkowy dt z warunku CFL
		dt = ( dx*dx*dy*dy / (2*vis * (dx*dx + dy*dy)) ) * 0.2;

		gifTime = 0.0;

		timePassed = 0.0;

		frames = 0;

		dtER = false;

		lidSpeed = 0.0;

		sU = false;
		sV = false;
		sP = false;
		sD = false;
};

System::~System(){
	cout<< "System: Destruktor"<< endl;
}

void System::addBorder(bool slip) {
	Type borderType;
	cout << nx << " " << ny << endl;
	if(slip)
	{
		borderType = border;
	}
	else {
		borderType = border_noslip;
	}

	for(int i=0; i<nx; ++i)
	{
		F1[i][0] = borderType; 
		F1[i][ny-1] = borderType; 

	}

	for(int j=1; j<ny-1; ++j)
	{
		F1[0][j] = borderType;
		F1[nx-1][j] = borderType;
	}
};

void System::addParticles(int n){
	double gap = 1.0/n; // przedział miedzy cząstkami
	double start = gap/2.0; // położenie pierwszej cząstki w komórce

	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){
			if(F1[i][j].fluid()){
				for(int k=0; k<n; k++)
				{
					for(int l=0; l<n; l++)
					C1.push_back(Particle( (i + start + k*gap) * dx, (j + start + l*gap) * dy, 0.0, 0.0, F1[i][j].getColor()) );
				}

			}
		}
	}
}

void System::fluidWall(int left, int right, int down, int up, int color){
	for(int i=left; i<right; ++i){
		for(int j=down; j<up; ++j){
			F1[i][j] = full;
			F1[i][j].setColor(color);
		}
	}
}


void System::solidWall(int left, int right, int down, int up, bool slip){ 
	Type borderType;
	if(slip)
	{
		borderType = border;
	}
	else {
		borderType = border_noslip;
	}

	for(int i=left; i<right; ++i){
		for(int j=down; j<up; ++j){
			F1[i][j] = borderType;
		}
	}
}

void System::fluidCircle(int x, int y, int r, int color){
	for(int i=0; i<nx; ++i){
		for(int j=0; j<ny; ++j){
			if (pow(i-x, 2) + pow(j-y, 2) < pow(r, 2)) 
			{
				F1[i][j] = full;
				F1[i][j].setColor(color);
			}
		}
	}
}

void System::addOut(int pos)
{
	F1[pos][0] = out;
}

void System::setLidSpeed(double s){
	lidSpeed = s;
}

void System::fluidSlope(int left, int right, int down, int up, int n, int color){
	double gap = 1.0/n; // przedział miedzy cząstkami
	double start = gap/2.0; // położenie pierwszej cząstki w komórce

	double ws = (right-left)*dx;
	double hs = (up-down)*dy;
	for(int i=left; i<right; ++i){
		for(int j=down; j<up; ++j){
			for(int k=0; k<n; k++){
				for(int l=0; l<n; l++){
					double xpos = (i + start + k*gap) * dx;
					double ypos = (j + start + l*gap) * dy;
					double d1 = std::sqrt(std::pow((xpos - left*dx)/ws, 2) + std::pow((ypos - down*dy)/hs, 2));
					double d2 = std::sqrt(std::pow((xpos - right*dx)/ws, 2) + std::pow((ypos - up*dy)/hs, 2));
						
					if(d1 < d2)
					{
						C1.push_back(Particle( xpos, ypos, 0.0, 0.0, color));

					}

				}
			}
		}
	}
}

void System::bubble(int x, int y, int r, int color){
	for(int i=0; i<nx; ++i){
		for(int j=0; j<ny; ++j){
			if (pow(i-x, 2) + pow(j-y, 2) > pow(r, 2)) 
			{
				if(!F1[i][j].boundary())
				{
					F1[i][j] = full;
					F1[i][j].setColor(color);
				}

			}
		}
	}
}


void System::saveF(){
  	ofstream myfile;
  	myfile.open ("../datafiles/komorki.dat");
	
	for (int i = 0; i<nx; ++i){
		for(int j=0; j<ny; ++j){
			myfile << i << " " << j << " " << F1[i][j].get()<< endl; 
		}
	}
  myfile.close();

};


void System::saveC() const{ 
  	ofstream myfile;
	if(frames==0) // zaczynamy nowy plik
  		myfile.open ("../datafiles/particles.dat");
	else // dopisujemy
		myfile.open ("../datafiles/particles.dat", std::ios_base::app); 
	
	for(const Particle& i : C1) 
		myfile << frames << " " << i.printToFile() << endl;
	myfile << endl;
	myfile << endl;
 	myfile.close();
}

void System::saveMatrix(Matrix<double> M, string filename)
{  	
	ofstream myfile;
	if(frames==0) // zaczynamy nowy plik
  		myfile.open ("../datafiles/" + filename);
	else // dopisujemy
		myfile.open ("../datafiles/" + filename, std::ios_base::app);
	
	for (int i = 0; i<nx; ++i){
		for(int j=0; j<ny; ++j){
			myfile << frames << " " << i << " " << j << " " << M[i][j]<< endl; 
		}
	}
	myfile << endl;
	myfile << endl;
 	myfile.close();
}

void System::saveU(){
	saveMatrix(U1, "velocitiesU.dat");
}

void System::saveV(){
	saveMatrix(V1, "velocitiesV.dat");
}

void System::saveP(){
	saveMatrix(P1, "pressures.dat");
}

void System::saveD(){
	saveMatrix(D1, "divergences.dat");
}

void System::saveMap()
{
	ofstream myfile;
  	myfile.open ("../datafiles/map.dat");
	
	for (int i = 0; i<nx; ++i){
		for(int j=0; j<ny; ++j){
			// if(F1.at(i, j).boundary() && F1.at(i, j)!=out)
			if(F1[i][j].boundary() && F1[i][j]!=out)
				myfile << (double(i)/double(nx)) * width << " " << (double(j)/double(ny)) * height << " " <<(double(i)/double(nx)) * width + dx<< " " << (double(j)/double(ny)) * height + dy << endl;
		}
	}
  myfile.close();
}

void System::print() const
{
    cout << "System: dt startowe: " << dt << endl;

}

void System::iteration(double dt, Matrix<double> &U, Matrix<double> &V, Matrix<double> &P,
	Matrix<double> &D, Matrix<Cell> &F, Matrix<double> &Z, std::vector<Particle> &C){

	//// --- SET UP ---

	// tablice na predkosci w kroku czasowym n+1
	Matrix<double> Un1;
	Matrix<double> Vn1;

	Un1.init(nx, ny);
	Vn1.init(nx, ny);

	//// --- PRZEFLAGOWYWANIE KOMOREK ---

	// sumy predkości cząstek w poszczególnych komórkach
	Matrix <double> U_c, V_c; 
	U_c.init(nx, ny);
	V_c.init(nx, ny);

 	// tablica pomocnicza mowiaca czy w danej komorce sa jakiekolwiek czastki i ile ich jest
	Matrix<int> Is_any_c; 
	Is_any_c.init(nx, ny); // defaltowo 0 -> false, brak czastek

	
	// zliczanie czastek i ich prędkości
	for(Particle& c : C)
	{
		int ix = floor(c.getX() / dx);
		int iy = floor(c.getY() / dy);

		Is_any_c[ix][iy] += 1;
		U_c[ix][iy]+= c.getUk();
		V_c[ix][iy]+= c.getVk();
	}

	// pierwsza poprawka do flagowania = wypelnianie 'od nowa'
	for(int i=1; i<nx-1; i++){  
		for(int j=1; j<ny-1; j++){
			if(Is_any_c[i][j] >0){ 
				// usrednianie predkosci - wykonujemy tutaj, później nie trzeba
				if(Is_any_c[i][j] > 0)
					U_c[i][j] /= double(Is_any_c[i][j]);
				else
					U_c[i][j] = 0.0;

				if(V_c[i][j] > 0)
					V_c[i][j] /= double(Is_any_c[i][j]);
				else
					V_c[i][j] = 0.0;

				if(F[i][j] == emptyc) // jest cząstka a było pusto: pusta -> powierzchniowa
				{
					F[i][j] = surface;
					if(F[i+1][j] == emptyc)
						U[i][j] = U_c[i][j];
					if(F[i-1][j] == emptyc)
						U[i-1][j] = U_c[i][j];
					if(F[i][j+1] == emptyc)
						V[i][j] = V_c[i][j];
					if(F[i][j-1] == emptyc)
						V[i][j-1] = V_c[i][j];
				}
				else{ // jest cząstka a i tak była: dajemy pełna, w następnej pętli powierzchniowe będą obsłużone
					F[i][j] = full;
				}

			}
			else{ // gdy nie ma cieczy, ciśnienie jest 0!
				P[i][j] =  0.0;
			}
		}
	}
	

	// druga poprawka do flagowania - ustawianie pelna -> powierzchniowa
	for(int i=1; i<nx; i++){
		for(int j=1; j<ny; j++){
			if(F[i][j] == full){
				if(F[i-1][j] == emptyc || F[i+1][j] == emptyc ||
					F[i][j-1] == emptyc || F[i][j+1] == emptyc){
						F[i][j] =  surface;
				}
			}
		}
	} 

	//przeflagowywanie - powierzchniowa <-> pusta
	for(int i=0; i<nx; ++i)
	{
		for(int j=0; j<ny; ++j)
		{
			if(!Is_any_c[i][j]) // nie ma czastki...
			{
				if(F[i][j] == surface) // ... a byla powierzchnia
				{
					F[i][j] = emptyc;
					P[i][j] = 0.0; // ważne by zerować ciśnienie w pustych!
					// zerowanie prędkości również pustych sąsiadów
					if(F[i+1][j]==emptyc) 
						U[i][j] = 0.0;
					if(F[i-1][j]==emptyc)
						U[i-1][j] = 0.0;
					if(F[i][j+1]==emptyc)
						V[i][j] = 0.0;
					if(F[i][j-1]==emptyc)
						V[i][j-1] = 0.0;
				}
			}
			else // jest czastka ...
			{
				// usrednianie predkosci - zrobione w pętli wcześniej

				if(F[i][j] == emptyc) // ... a było pusto: pusta -> powierzchniowa
				{
					F[i][j] = surface;
					if(F[i+1][j] == emptyc || F[i+1][j] == out)
						U[i][j] = U_c[i][j];
					if(F[i-1][j] == emptyc || F[i-1][j] == out)
						U[i-1][j] = U_c[i][j];
					if(F[i][j+1] == emptyc || F[i][j+1] == out)
						V[i][j] = V_c[i][j];
					if(F[i][j-1] == emptyc || F[i][j-1] == out)
						V[i][j-1] = V_c[i][j];

					// ciśnienie dla powierzchniowej zostaje 0 jak w pustej
				}

			}
		}
	}

	//flagowanie pelna <-> powierzchniowa
	for(int i=0; i<nx; ++i)
	{
		for(int j=0; j<ny; ++j)
		{
			if(Is_any_c[i][j])
			{
				if(F[i][j] == full)
				{
					// ktorykolwiek sasiad pusty : pelna -> powierzchniowa
					if(F[i+1][j]==emptyc || F[i-1][j] == emptyc || F[i][j+1] == emptyc || F[i][j-1] == emptyc)
					{
						F[i][j] = surface;
						P[i][j] = 0.0; // powierzchniowa: cisnienie 0
					}
				}
			 	if (F[i][j] == surface) 
				{
					// sprawdzamy czy wszyscy sąsiedzi są pełni

					// kwestia do sprawdzenia: w pracy Harlova spradzamy czy każdy sąsiad jest pełny, i tylko wtedy zmieniamy pow -> pełna
					// teoretycznie możliwe jest żeby któryś sąsiad był powierzchniowy, wtedy powinniśmy porównywać z F[][].fluid()
					// możliwe jest że jest to celowy zabieg by poprawić stabilność w programie (zbyt 'pochopne' przeflagowywanie)

					double p_tmp = 0.0;
					int full_count = 0; // full counter or sth
					
					if(F[i+1][j] == full){ 
						p_tmp += P[i+1][j];
						full_count++;
					}
					if(F[i-1][j] == full)
					{
						p_tmp += P[i-1][j];
						full_count++;
					}
					if(F[i][j+1] == full)
					{
						p_tmp += P[i][j+1];
						full_count++;
					}
					if(F[i][j-1] == full)
					{
						p_tmp += P[i][j-1];
						full_count++;
					}

						
					// wszyscy sasiedzi pelni: pow -> pelna
					// if( Is_any_c[i+1][j] && Is_any_c[i-1][j] && Is_any_c[i][j+1] && Is_any_c[i][j-1]) // ten warunek powinien być równoważny z count==4
					if(full_count==4)
					{
						F[i][j] = full;
						P[i][j] = p_tmp/4.0;
					}
						
				}
			}
		}
	}

	
	//// Usuwanie cząstek na wyjściu

	std::vector<Particle>::size_type size = C.size();

    for (std::vector<Particle>::size_type i = 0; i < size; ++i)
    {
		Particle& c = C[i]; 
		int ix = floor(c.getX() / dx);
		int iy = floor(c.getY() / dy); 

		if(F[ix][iy].boundary()) // == out // usuwamy każdą cząstkę która trafi na out
		{
 			C.erase(C.begin() + i);
		}
	}

	//// Ustawienie predkosci na wieczku
	for(int i=0; i<nx; i++)
		U[i][ny-1] = lidSpeed;
	
	//// --- DYWERGENCJA ---

	for(int i=0; i<nx; ++i)
	{
		for(int j=0; j<ny; ++j)
		{
			if(F[i][j]==full) 
			{
				D[i][j] = ( U[i][j] - U[i-1][j])/dx + (V[i][j] - V[i][j-1] )/dy;
			}
			else if(F[i][j] == surface)
			{
				D[i][j] = 0.0;
			}

		}
	}

	//// --- ROZKŁAD ŹRÓDEŁ --- 
	
	double sumZ = 0.0; // suma źródeł, używana później do warunku zbieżności ciśnień
	for(int j=1; j<ny-1; ++j)
	{
		for(int i=1; i<nx-1; ++i)
		{
			if(F[i][j] == full)
			{
				double S1 = 0.0;
				double S2 = 0.0;

				// pierwsza czesc : obsługa narożników
				if(COR(F, i+1, j))
				{
					if(OB(F, i+1, j+1)){
						S2 -= U[i][j+1]/dx;
					}
					else{
						S2 -= U[i][j-1]/dx;
					}
				}
				if(COR(F, i-1, j))
				{
					if(OB(F, i-1, j+1)){
						S2 += U[i-1][j+1]/dx;
					}
					else{
						S2 += U[i-1][j-1]/dx;
					}
				}
				if(COR(F, i, j+1))
				{
					if(OB(F, i+1, j+1)){
						S1 -= V[i+1][j]/dy;
					}
					else{
						S1 -= V[i-1][j]/dy;
					}
				}
				if(COR(F, i, j-1))
				{
					if(OB(F, i+1, j-1)){
						S1 += V[i+1][j-1]/dy;
					}
					else{
						S1 -= V[i-1][j-1]/dy;
					}
				}

				// druga czesc
				if(!OB(F, i, j))
				{
					if(COR(F, i+1, j-1))
					{
						if(F[i+1][j-1].slip())
						{
							S2 += U[i][j]/dx;
							S1 -= V[i][j-1]/dy;
						}
						else{
							S2 -= U[i][j]/dx;
							S1 += V[i][j-1]/dy;
						}
					}
					if(COR(F, i-1, j-1))
					{
						if(F[i-1][j-1].slip())
						{
							S2 -= U[i-1][j]/dx;
							S1 -= V[i][j-1]/dy;
						}
						else{
							S2 += U[i-1][j]/dx;
							S1 += V[i][j-1]/dy;
						}
					}
					if(COR(F, i+1, j+1))
					{
						if(F[i+1][j+1].slip())
						{
							S2 += U[i][j]/dx;
							S1 += V[i][j]/dy;
						}
						else{
							S2 -= U[i][j]/dx;
							S1 -= V[i][j]/dy;
						}
					}
					if(COR(F, i-1, j+1))
					{
						if(F[i-1][j+1].slip())
						{
							S2 -= U[i-1][j]/dx;
							S1 += V[i][j]/dy;
						}
						else{
							S2 += U[i-1][j]/dx;
							S1 -= V[i][j]/dy;
						}
					}
				}
				/////////////////////

				double Z1 = 0.0;
				double Z4 = 0.0;

				// Z1 oraz Z4
				if(!F[i+1][j].boundary()) // po prawej nie ma brzegu
				{
					Z1+= (pow(( U[i][j] + U[i+1][j] ) / 2, 2) - pow((U[i-1][j] + U[i][j])/2, 2) ) / pow(dx, 2);
					Z4+= (D[i+1][j] - D[i][j]) / pow(dx, 2);
				}
				if(!F[i-1][j].boundary()) // po lewej nie ma brzegu
				{
					Z1+= (pow(( U[i-1][j] + U[i-2][j] ) / 2, 2) - pow( (U[i-1][j] + U[i][j])/2, 2) ) / pow(dx, 2);
					Z4+= (D[i-1][j] - D[i][j]) / pow(dx, 2);
				}
				if(!F[i][j+1].boundary()) // gorny nie jest brzegiem
				{
					Z1+= (pow(( V[i][j] + V[i][j+1] ) / 2, 2) - pow((V[i][j-1] + V[i][j])/2, 2) ) / pow(dy, 2);
					Z4+= (D[i][j+1] - D[i][j]) / pow(dy, 2);
				}
				if(!F[i][j-1].boundary()) // dolny nie jest brzegiem
				{
					Z1+= (pow(( V[i][j-1] + V[i][j-2] ) / 2, 2) - pow((V[i][j-1] + V[i][j])/2, 2) ) / pow(dy, 2);
					Z4+= (D[i][j-1] - D[i][j]) / pow(dy, 2);
				}

				// Z2 został zastąpiony przez Z5 poniżej - zostawiam do wglądu, brakowało tu warunków
				// Z2
				// double Z2 = (U[i][j] + U[i][j+1]) * (V[i][j] + V[i+1][j]);
				// Z2 += (U[i-1][j] + U[i-1][j-1]) * (V[i][j-1] + V[i-1][j-1]);
				// Z2 -= (U[i][j] + U[i][j-1]) * (V[i][j-1] + V[i+1][j-1]);
				// Z2 -= (U[i-1][j] + U[i-1][j+1]) * (V[i-1][j] + V[i][j]);
				// Z2 /= (2 * dx * dy);

				// Z3
				double Z3 = -D[i][j]/dt; 

				// Z5
				double Z5 = 0.0;
				if(!URON(F, i, j)){
					Z5 += (U[i][j] + U[i][j+1]) * (V[i][j] + V[i+1][j]);
				}
				if(!URON(F, i-1, j-1))
				{
					Z5+= (U[i-1][j] + U[i-1][j-1]) * (V[i][j-1] + V[i-1][j-1]);
				}
				if(!URON(F, i, j-1))
				{
					Z5-= (U[i][j] + U[i][j-1]) * (V[i][j-1] + V[i+1][j-1]);
				}
				if(!URON(F, i-1, j))
				{
					Z5-= (U[i-1][j] + U[i-1][j+1]) * (V[i-1][j] + V[i][j]);
				}
				Z5 /= (2 * dx * dy);

				// suma
				Z[i][j] = Z1 + Z3 + Z5 + S1/pow(dx, 2) + S2/pow(dy, 2) - vis*Z4; 
				sumZ += fabs(Z[i][j]); // zwiększenie sumy źródeł
			}
			else{ 
				// ewentualna obsługa źródeł dla niepełnej komórki
			}

		}
	}

	//// --- RELAKSACJA CISNIEN ---

	DiffPMax = 2.0 * EpsP;
	double sumP = 2.0 * EpsP; // suma ciśnień do drugiego warunku
	while (DiffPMax > EpsP && sumP > EpsP)
	{
		DiffPMax = 0.0;
		sumP = 0.0;
		for(int i=1; i<nx-1; ++i)
		{
			for(int j=1; j<ny-1; ++j)
			{
				if(F[i][j].fluid()) // wg algorytmu tylko pełna, ale z tym działa lepiej (?)
				{
					//ustal warunki brzegowe na P
					if(OB(F, i, j))
						boundaryConditionsP(i, j, F, P, U, V);
					
					double oldP = P[i][j];
					// nowe P
					P[i][j] = ( P[i+1][j] +P[i-1][j]) * pow(dy, 2) + ( P[i][j+1] + P[i][j-1]) * pow(dx, 2);
					P[i][j] += Z[i][j] * pow(dx, 2) * pow(dy, 2);
					P[i][j] /= (2 * (pow(dx, 2) + pow(dy, 2))); 

					// roznica ciśnien
					double DiffP = fabs(P[i][j] - oldP);

					// normalizacja roznicy cisnien
					DiffP /= (fabs(P[i][j]) + fabs(oldP) + pow((U[i][j] + U[i-1][j])/2.0, 2) +
								pow((V[i][j] + V[i][j-1])/2.0, 2) + fabs(gy)*width + fabs(gx)*width); 
	
					DiffPMax = DiffPMax > DiffP? DiffPMax : DiffP; 

					sumP += fabs(P[i+1][j] + P[i-1][j] - 2*P[i][j])/pow(dx, 2);
					sumP += fabs(P[i][j+1] + P[i][j-1] - 2*P[i][j])/pow(dy, 2);
				}

			}
		}
		// obliczenie sumy ciśnień po całych 2 pętlach
		sumP = fabs(sumP/sumZ -1.0);
	}

	//// --- PRĘDKOŚCI ---
	for(int i=0; i< nx; i++)
	{
		for(int j=0; j<ny; j++)
		{
			if(F[i][j].fluid())
			{
				// część dla U (#1)
				if(F[i+1][j]==emptyc)
				{
					Un1[i][j] = U[i][j] + gx*dt;
				}
				else if(F[i+1][j].boundary())
				{
					Un1[i][j] = U[i][j];
				}
				else 
				{
					double ui12jp1 = U[i][j+1];
					double ui12jm1 = U[i][j-1];

					// używanie sąsiada w razie gdy pusty
					if(F[i][j+1] == emptyc  && F[i+1][j+1] == emptyc)
					{
						ui12jp1 = U[i][j];
					}
					if(F[i][j-1] == emptyc  && F[i+1][j-1] == emptyc)
					{
						ui12jm1 = U[i][j];
					}

					double SU1 =( U[i+1][j] + U[i-1][j] -2*U[i][j])/pow(dx, 2); 
					SU1 += (ui12jp1 + ui12jm1 - 2*U[i][j])/pow(dy, 2);
					SU1 *= vis;
					SU1 += (P[i][j] - P[i+1][j])/dx + gx;
					
					double SU2 = 0.0;
					if(!URON(F, i, j))
						SU2 -= ((U[i][j] + ui12jp1)/2.0) * ((V[i][j] + V[i+1][j])/2.0);
					if(!URON(F, i, j-1))
						SU2 += ((ui12jm1 + U[i][j])/2.0) * ((V[i][j-1] + V[i+1][j-1])/2.0);


					double SU3  = pow((U[i-1][j] + U[i][j])/2.0, 2);
					SU3 -= pow(( U[i+1][j] + U[i][j] )/2.0, 2);

					Un1[i][j] = U[i][j] + dt*(SU1 + SU2/dy + SU3/dx);
				}

				// część dla V 
				if(F[i][j+1]==emptyc)
				{
					Vn1[i][j] = V[i][j] + gy*dt;
				}
				else if(F[i][j+1].boundary())
				{
					Vn1[i][j] = V[i][j];
				}
				else 
				{
					double vip1j12 = V[i+1][j];
					double vim1j12 = V[i-1][j];
					// uzywanie sąsiada w razie gdy pusty
					if(F[i+1][j] == emptyc  && F[i+1][j+1] == emptyc)
					{
						vip1j12 = V[i][j];
					}
					if(F[i-1][j] == emptyc  && F[i-1][j+1] == emptyc)
					{
						vim1j12 = V[i][j];
					}

					double SV1 = (vip1j12 + vim1j12 - 2*V[i][j]) /pow(dx, 2); 
					SV1 += (V[i][j+1] + V[i][j-1] - 2*V[i][j]) /pow(dy, 2);
					SV1 *= vis;
					SV1 += (P[i][j] - P[i][j+1])/dy + gy;
				
					double SV2 = 0.0;
					if(!URON(F, i, j))
						SV2 -= ( (U[i][j] + U[i][j+1])/ 2.0) * ((V[i][j] + vip1j12) / 2.0);
					if(!URON(F, i, j-1))
						SV2 += ((U[i-1][j] + U[i-1][j+1])/2.0) * ((vim1j12 + V[i][j])/2.0);


					double SV3 = pow((V[i][j-1] + V[i][j])/2.0, 2);
					SV3 -= pow((V[i][j+1] + V[i][j])/2.0, 2);

					Vn1[i][j] = V[i][j] + dt*(SV1 + SV2/dx + SV3/dy);

				}
			}
		}
	}

	for(int i=0; i<nx; i++)
	{
		for(int j=0; j<ny; j++)
		{
			if(F[i][j]==surface) 
			{
				if(F[i-1][j] == emptyc)
					Un1[i-1][j] = Un1[i][j];
				else if(F[i+1][j] == emptyc)
					Un1[i][j] = Un1[i-1][j];

				if(F[i][j-1] == emptyc)
					Vn1[i][j-1] = Vn1[i][j];
				else if(F[i][j+1] == emptyc)
					Vn1[i][j] = Vn1[i][j-1];
			}
			else{
				if(F[i][j].boundary())
				{
					if(OB(F, i+1, j))
						Un1[i][j] = U[i][j];
					if(OB(F, i, j+1))
						Vn1[i][j] = V[i][j];
				}
			}
		}
	}

	
	for(int i=0; i<nx; i++)
	{
		for(int j=0; j<ny; j++)
		{
			if(OB(F, i, j))
			{
				if(!(F[i][j] == emptyc))
				{
					// warunki na OUT ...
					
					if(F[i+1][j].boundary()){
						if(F[i+1][j]==out){
							Un1[i][j] = Un1[i-1][j] - (dx/dy) * (Vn1[i][j] - Vn1[i][j-1]);
						} else {
							Un1[i][j] = U[i][j];
						}
					}
					if(F[i-1][j].boundary()){
						if(F[i-1][j]==out){
							Un1[i-1][j] = Un1[i][j] + (dx/dy) * (Vn1[i][j] - Vn1[i][j-1]);
						} else {
							Un1[i-1][j] = U[i-1][j];
						}
					}
					if(F[i][j+1].boundary()){
						if(F[i][j+1]==out){
							Vn1[i][j] = Vn1[i][j-1] - (dy/dx) * (Un1[i][j] - Un1[i-1][j]);
						} else {
							Vn1[i][j] = V[i][j];
						}
					}
					if(F[i][j-1].boundary()){
						if(F[i][j-1]==out){
							Vn1[i][j-1] = Vn1[i][j] + (dy/dx) * (Un1[i][j] - Un1[i-1][j]);
						} else {
							Vn1[i][j-1] = V[i][j-1];
						}	
					}
				}
			}
		}
	}

	// warunki brzegowe dla prędkości
	for(int i=0; i<nx; i++)
	{
		for(int j=0; j<ny; j++)
		{
				// scianki ograniczajace
				if(!F[i][j].boundary())
				{
					if (F[i+1][j].boundary()) // scianka po prawej
					{
						if(F[i+1][j].slip())
						{
							Vn1[i+1][j] = Vn1[i][j];
							Un1[i][j] = 0.0;
						}
						else if(F[i+1][j].noSlip())
						{
							Vn1[i+1][j] = -1.0 * Vn1[i][j];
							Un1[i][j] = 0.0;
						}
						else if(F[i+1][j]==out)
						{
							Vn1[i+1][j] = Vn1[i][j];
						}
					}
					
					if(F[i-1][j].boundary()) // scianka po lewej
					{
						if(F[i-1][j].slip())
						{
							Vn1[i-1][j] = Vn1[i][j];
							Un1[i-1][j] = 0.0;
						}
						else if(F[i-1][j].noSlip())
						{
							Vn1[i-1][j] = -1.0 * Vn1[i][j];
							Un1[i-1][j] = 0.0;
						}
						else if(F[i-1][j] == out)
						{
							Vn1[i-1][j] = Vn1[i][j];
						}
					}
					if(F[i][j-1].boundary()) // scianka na dole
					{
						if(F[i][j-1].slip())
						{
							Vn1[i][j-1] = 0.0;
							Un1[i][j-1] = Un1[i][j];
						}
						else if(F[i][j-1].noSlip())
						{
							Vn1[i][j-1] = 0.0;
							Un1[i][j-1] = -1.0 * Un1[i][j];
						}
						else if(F[i][j-1]==out)
						{
							Un1[i][j-1] = Un1[i][j];
						}
					}
					if(F[i][j+1].boundary()) // scianka na gorze
					{
						if(F[i][j+1].slip()) 
						{
							Vn1[i][j] = 0.0; 
							Un1[i][j+1] = Un1[i][j];
						}
						else if(F[i][j+1].noSlip())
						{
							Vn1[i][j] = 0.0; 
							Un1[i][j+1] = -1.0 * Un1[i][j];
						}
						else if(F[i][j+1]==out)
						{
							Un1[i][j+1] = Un1[i][j];
						}
					}
				}
		}
	}


	//// --- warunki dla powierzchni cieczy  (ciecz - atmosfera) ---
	for(int i=0; i<nx; ++i)
	{
		for(int j=0; j<ny; ++j)
		{
			if(F[i][j]==surface)
			{
				emptyNeighbors(i, j, F);
				if(neighSum==1) //jedna pusta
				{
					if(neighbors[0]){ // north
						Vn1[i][j] = Vn1[i][j-1] + (dy/dx) * (Un1[i-1][j] - Un1[i][j]);
					}else if(neighbors[1]){ //east
						Un1[i][j] = Un1[i-1][j] + (dx/dy) * (Vn1[i][j-1] - Vn1[i][j]);
					}else if(neighbors[2]){ //south
						Vn1[i][j-1] = Vn1[i][j] + (dy/dx) * (Un1[i][j] - Un1[i-1][j]);
					}else if(neighbors[3]){//west
						Un1[i-1][j] = Un1[i][j] + (dx/dy) * (Vn1[i][j] - Vn1[i][j-1]);
					}
				}
				else if(neighSum==2)
				{
					if(neighbors[0]!=neighbors[2]) // rozne naprzeciw - sasiedzi
					{
						if(neighbors[0]){
							Vn1[i][j] = Vn1[i][j-1];
						}else if(neighbors[2]){
							Vn1[i][j-1] = Vn1[i][j];
						}

						if(neighbors[1]){
							Un1[i][j] = Un1[i-1][j];
						}else if(neighbors[3]){
							Un1[i-1][j] = Un1[i][j];
						}
					}
					else // naprzeciw siebie puste
					{
						if(neighbors[0]) // puste gora-dol
						{
							Un1[i-1][j] += dt*gx;
							Un1[i][j] += dt*gx;
						}
						else if(neighbors[1]) // puste lewy-prawy
						{
							Vn1[i][j] += dt*gy;
							Vn1[i][j-1] += dt*gy;
						}
					}
				}
				else if(neighSum==3)
				{
					if(!neighbors[0])
					{
						Un1[i-1][j] += dt*gx;
						Un1[i][j] += dt*gx;
						Vn1[i][j-1] = Vn1[i][j];
					}
					else if(!neighbors[1])
					{
						Vn1[i][j] += dt*gy;
						Vn1[i][j-1] += dt*gy;
						Un1[i-1][j] = Un1[i][j];
					}
					else if(!neighbors[2])
					{
						Un1[i-1][j] += dt*gx;
						Un1[i][j] += dt*gx;
						Vn1[i][j] = Vn1[i][j-1];
					}
					else if(!neighbors[3])
					{
						Vn1[i][j] += dt*gy;
						Vn1[i][j-1] +=dt*gy;
						Un1[i][j] = Un1[i-1][j]; 
					}
				}
				else if(neighSum==4)
				{
					Un1[i-1][j] +=dt*gx;
					Un1[i][j] += dt*gx;
					Vn1[i][j] += dt*gy;
					Vn1[i][j-1] += dt*gy;
				}
			}
		}
	}
	// przepisanie nowych prędkości - robimy na samym końcu (wbrew pracy)
	U = Un1;
	V = Vn1;


	//// --- PRZESUNIECIE CZASTEK ZNACZONYCH --- 
	for(Particle& c : C)
	{
		double x = c.getX();
		double y = c.getY();

		int i = floor(x / dx);
		int j = floor(y / dy);

		// U
		double dpos = (y/dy) - j;

		int jp = dpos<0.5 ? j-1 : j;

		double Px = i - x/dx + 0.5;
		double Py = jp - y/dy + 1.0;

		double w1 = (0.5 + Px)*(0.5 - Py);
		double w2 = (0.5 - Px)*(0.5 - Py);
		double w3 = (0.5 + Px)*(0.5 + Py);
		double w4 = (0.5 - Px)*(0.5 + Py);

		double u1 = 0.0;
		double u2 = 0.0;
		double u3 = 0.0;
		double u4 = 0.0;

		if(fabs(U[i-1][jp+1])<=zeroEps){
			u1 = U[i-1][j];
		}else if(!COR(F, i-1, j)){
			u1 = U[i-1][jp+1];
		}

		if(fabs(U[i][jp+1]) <= zeroEps){
			u2 = U[i][j];
		}else if(!COR(F, i+1, j)){
			u2 = U[i][jp+1];
		}

		if(fabs(U[i-1][jp]) <= zeroEps){
			u3 = U[i-1][j];
		}else if (!COR(F, i-1, j)){
			u3 = U[i-1][jp];
		}

		if(fabs(U[i][jp])<= zeroEps){
			u4 = U[i][j];
		}else if(!COR(F, i+1, j)){
			u4 = U[i][jp];
		}

		double uk = w1*u1 + w2*u2 + w3*u3 + w4*u4;
		c.setUk(uk);
		double xn = x + uk*dt;
		c.setX(xn);
		//V
		dpos = (x/dx) - i;

		int ip = dpos<0.5 ? i-1 : i;

		Px = ip - x/dx + 1.0;
		Py = j - y/dy + 0.5;

		w1 = (0.5 + Px) * (0.5 - Py);
		w2 = (0.5 - Px) * (0.5 - Py);
		w3 = (0.5 + Px) * (0.5 + Py);
		w4 = (0.5 - Px) * (0.5 + Py);

		double v1 = 0.0;
		double v2 = 0.0;
		double v3 = 0.0;
		double v4 = 0.0;


		if(fabs(V[ip][j])<=zeroEps){
			v1 = V[i][j];
		}else if(!COR(F, i, j+1)){
			v1 = V[ip][j]; 
		} 

		if(fabs(V[ip+1][j]) <= zeroEps){ 
			v2 = V[i][j];
		}else if(!COR(F, i, j+1)){
			v2 = V[ip+1][j];
		}

		if(fabs(V[ip][j-1]) <= zeroEps){
			v3 = V[i][j-1];
		}else if(!COR(F, i, j-1)){
			v3 = V[ip][j-1];
		}

		if(fabs(V[ip+1][j-1])<= zeroEps){
			v4 = V[i][j-1];
		}else if(!COR(F, i, j-1)){
			v4 = V[ip+1][j-1];
		}

		double vk = w1*v1 + w2*v2 + w3*v3 + w4*v4;
		c.setVk(vk);
		double yn = y + vk*dt;
		c.setY(yn);

	}
	// koniec funkcji iteration
}

void System::boundaryConditionsP(int i, int j, Matrix<Cell> &F, Matrix<double> &P,
	Matrix<double> &U, Matrix<double> &V)
{
	// w razie rozbudowy do warunkow z noSlip dodać lub in, a w środku warunek na incor
	if(F[i][j+1].noSlip()) 
	{
		P[i][j+1] = P[i][j] + gy*dy + 2*vis*( V[i][j-1] - V[i][j] )/dy;
	}
	else if(F[i][j+1].slip())
	{
		P[i][j+1] = P[i][j] + gy*dy;
	}
	else if(F[i][j+1]==out)
	{
		double a = ((U[i-1][j] + U[i-1][j+1])/2) * ((V[i-1][j] + V[i][j])/2);
		double b = ((U[i][j] + U[i][j+1])/2) * ((V[i][j] + V[i+1][j])/2);

		P[i][j+1] = P[i][j] + (dy/dx) * (a - b);
	}

	if(F[i][j-1].noSlip())
	{
		P[i][j-1] = P[i][j] - gy*dy - 2*vis*( V[i][j] - V[i][j-1] )/dy;
	}
	else if(F[i][j-1].slip())
	{
		P[i][j-1] = P[i][j] - gy*dy;
	}
	else if(F[i][j-1]==out)
	{
		double a = ((U[i-1][j-1] + U[i-1][j])/2) * ((V[i-1][j-1] + V[i][j-1])/2);
		double b = ((U[i][j] + U[i][j-1])/2) * ((V[i][j-1] + V[i+1][j-1])/2);

		P[i][j-1] = P[i][j] - (dy/dx) * (a - b);
	}
	

	if(F[i+1][j].noSlip())
	{
		P[i+1][j] = P[i][j] + gx*dx + 2*vis*( U[i-1][j] - U[i][j] )/dx;
	}
	else if(F[i+1][j].slip())
	{
		P[i+1][j] = P[i][j] + gx*dx;
	}
	else if(F[i+1][j]==out)
	{
		double a = ((U[i][j-1] + U[i][j])/2) * ((V[i][j-1] + V[i+1][j-1])/2);
		double b = ((U[i][j] + U[i][j+1])/2) * ((V[i][j] + V[i+1][j])/2);

		P[i+1][j] = P[i][j] + (dx/dy) * (a - b);
	}

	if(F[i-1][j].noSlip())
	{
		P[i-1][j] = P[i][j] - gx*dx - 2*vis*( U[i][j] - U[i-1][j] )/dx;
	}
	else if(F[i-1][j].slip())
	{
		P[i-1][j] = P[i][j] - gx*dx;
	}
	else if(F[i-1][j]==out)
	{
		double a = ((U[i-1][j-1] + U[i-1][j])/2) * ((V[i][j-1] + V[i-1][j-1])/2);
		double b = ((U[i-1][j] + U[i-1][j+1])/2) * ((V[i][j] + V[i-1][j])/2);

		P[i-1][j] = P[i][j] - (dx/dy) * (a - b);
	}
}

void System::emptyNeighbors(int i, int j, Matrix<Cell> &F)
{
	neighbors[0] = (F[i][j+1] == emptyc); // gorna
	neighbors[1] = (F[i+1][j] == emptyc); // prawa
	neighbors[2] = (F[i][j-1] == emptyc); // dolna
	neighbors[3] = (F[i-1][j] == emptyc); // lewa

	neighSum = 0;
	for(int k=0; k<4; k++)
	{
		if(neighbors[k])
			neighSum++;
	}
}


bool System::fullNeighbors(Matrix<int> &Is_any_c, int i, int j)
{
	if(Is_any_c[i+1][j] == 0 || Is_any_c[i-1][j] == 0 || Is_any_c[i][j+1] == 0 || Is_any_c[i][j-1]==0)
		return false;
	return true; 
}

void System::step(){
	if(!dtER){
		stepCFL();
	}
	else {
		stepRE();
	}
}

void System::stepCFL(){
	iters++;
	// kolejna iteracja
	iteration(dt, U1, V1, P1, D1, F1, Z1, C1);

	timePassed += dt;

	if(timePassed >= gifTime) // przekroczony kolejny prog czasu do zapisu?
	{
		saveC();				// zapisz
		if(sU)
			saveU();
		if(sV)
			saveV();
		if(sP)
			saveP();
		if(sD)
			saveD(); 			
		gifTime += frameTime; 	// nowy próg 
		frames+=1; 				// zwiększ ilość klatek
	}

	// wyliczenie nowego dt
	newDtCFL();
}

void System::stepRE(){
	// primowane = glowne // Fp = F itd
	iters++;

	// kopie do iteracji raz całym dt
	Matrix<double> U2 = U1;
	Matrix<double> V2 = V1; 
	Matrix<double> P2 = P1;
	Matrix<double> D2 = D1;
	Matrix<Cell> F2 = F1;
	Matrix<double> Z2 = Z1;
	std::vector<Particle> C2 = C1;

	// zachowuje wartosci z poczatku iteracji zeby do nich wrocic jesli krok nie zaakceptowany
	Matrix<double> U1tmp = U1; 
	Matrix<double> V1tmp = V1;
	Matrix<double> P1tmp = P1;
	Matrix<double> D1tmp = D1;
	Matrix<Cell> F1tmp = F1;
	Matrix<double> Z1tmp = Z1;
	std::vector<Particle> C1tmp = C1;

	//2 razy 1/2 dt na wersjach głównych
	iteration(dt*0.5, U1, V1, P1, D1, F1, Z1, C1);
	iteration(dt*0.5, U1, V1, P1, D1, F1, Z1, C1);

	// 1 raz dt na wersjach kopiowanych
	iteration(dt, U2, V2, P2, D2, F2, Z2, C2); 

	// liczenie nowego dt
	bool accept = newDtRE(U2, V2);

	if(accept){ // jesli krok zaakceptowany, to wstawiamy wynik po 2 krokach (bo jest dokladniejszy)
		// nie trzeba kopiować, wszystko jest ok w wersjach z _1
		timePassed += dt;
		
		if(timePassed >= gifTime) // przekroczony kolejny prog czasu do zapisu?
		{
			saveC(); 				// zapisz
			if(sU)
				saveU();
			if(sV)
				saveV();
			if(sP)
				saveP();
			if(sD)
				saveD(); 
			gifTime += frameTime; 	// nowy próg 
			frames+=1; 				// zwiększ ilość klatek
		}
	}
	else{ // jesli krok nie zaakceptowany, to wracamy do poczatkowych
		U1 = U1tmp;
		V1 = V1tmp;
		P1 = P1tmp;
		D1 = D1tmp;
		F1 = F1tmp;
		Z1 = Z1tmp;
		C1 = C1tmp;
	}
}

bool System::newDtRE(Matrix<double> U2, Matrix<double> V2){
	double S = 0.8; 	// parametr < 1
	double T = 0.007; 	// tolerancja błędu - ma duzy wplyw na krok czasowy
	double p = 1.0; 	// rzad dokladnosci metody

	double dn = -1;
	double du, dv;

	// liczymy maksymalna różnicę miedzy predkosciami
	for(int i=0; i< nx; i++)
	{
		for(int j=0; j<ny; j++)
		{
			
			du = fabs(U1[i][j] - U2[i][j]); // 
			dv = fabs(V1[i][j] - V2[i][j]);
			if(du > dn)
				dn = du;
			if(dv > dn)
				dn = dv;
		}
	}

	// liczymy nowe dt  na bazie dn
	dt = pow(S * (T/dn), 1.0/(p+1.0)) * dt;

	// akceptujemy jedynie gdy dn jest mniejsze niz nasza zadana tolerancja bledu
	return T > dn;


}

void System::newDtCFL(){
	double u_max = 0.00;
	double v_max = 0.00;
	for(int i=0; i< nx; i++)
	{
		for(int j=0; j<ny; j++)
		{
			if(U1[i][j] > u_max)
				u_max = U1[i][j];
			if(V1[i][j] > v_max)
				v_max = V1[i][j];
		}
	}

	double lam_cfl = 0.5; // 0.6

	double new_dt = dx*dx*dy*dy / (2*vis * (dx*dx + dy*dy));

	if(u_max > 0.0 || v_max > 0.0)
	{
		double dt_u = lam_cfl * dx/u_max;

		double dt_v = lam_cfl * dy/v_max;
	
		if (dt_u < new_dt)
			new_dt = dt_u;
		if(dt_v < new_dt)
			new_dt = dt_v;

		// razy 2 dziel przez 2

		// parametr stabilizacji
		new_dt = 0.2 * new_dt; // 0.2
		if (dt <= new_dt) // moze byc wiekszy krok
		{
			// dt = (dt + new_dt)/2.0; 
			dt =  new_dt;
		}
		else{
			dt =  new_dt;
		}
	}
	else{
		dt = new_dt;
	}
	
}

bool System::URON(Matrix<Cell> F, int i, int j)
{
	if(i>0 && i < nx && j>0 && j < ny)
	{
		if(F[i+1][j+1].boundary())
			return true;
	}
	return false;
	
}

bool System::OB(Matrix<Cell> F, int i, int j)
{
	if(F[i][j].fluid() || F[i][j]==emptyc)
	{
		if(i>0 && i < nx && j>0 && j < ny)
		{
			if(F[i+1][j].boundary() || F[i-1][j].boundary() || F[i][j-1].boundary() || F[i][j+1].boundary())
				return true;
		}
	}
	return false;
}

bool System::COR(Matrix<Cell> F, int i, int j)
{
	if(i>0 && i < nx && j>0 && j < ny)
	{
		// brzeg _| ()
		if(F[i][j].boundary() && F[i-1][j].boundary() && F[i][j+1].boundary() && !F[i+1][j].boundary() && !F[i][j-1].boundary())
			return true;
		// brzeg |_
		if(F[i][j].boundary() && !F[i-1][j].boundary() && F[i][j+1].boundary() && F[i+1][j].boundary() && !F[i][j-1].boundary())
			return true;
		// brzeg -|
		if(F[i][j].boundary() && F[i-1][j].boundary() && !F[i][j+1].boundary() && !F[i+1][j].boundary() && F[i][j-1].boundary())
			return true;
		// brzeg |-
		if(F[i][j].boundary() && !F[i-1][j].boundary() && !F[i][j+1].boundary() && F[i+1][j].boundary() && F[i][j-1].boundary())
			return true;
	}
		return false;
}

void System::setER(){
	dtER = true;
}

void System::setCFL(){
	dtER = false;
}

double System::getTimePassed(){
	return timePassed;
}

void System::saveUVPD(bool u, bool v, bool p, bool d)
{
	sU = u;
	sV = v;
	sP = p;
	sD = d;
}

