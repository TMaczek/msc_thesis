#pragma once

#include "Particle.h"
#include "Cell.h"
#include "Matrix.h"
#include <vector>

/**
 * @brief 
 * Klasa reprezentująca system w którym znajduje się ciecz: kontener i tablice dla siatek Eulera dla wszystkich parametrów.
 */
class System{
    private:
		// wymiar poziomy kontenera
		double width;
		// wymiar pionowy kontenera
		double height;

		// ilosc komorek siatki w kierunku osi x 
		int nx;
		// ilosc komorek siatki w kierunku osi y
		int ny;

		// rozmiar kroku poziomego na siatce (= width/nx)
		double dx; 
		// rozmiar kroku pionowego na siatce (= height/ny)
		double dy; 

		// sily zewnetrzne w kierunku osi x
		double gx;
		// sily zewnetrzne w kierunku osi y (najczesciej: sila grawitacji)
		double gy; 

		// lepkosc kinamatyczna cieczy
		double vis; 

		// krok czasowy 
		double dt; 

		// ilosc minionego czasu symulacji
		double timePassed;

		// przedzial czasowy co ile zapisywać wyniki do pliku (żeby w animacji był stały)
		double frameTime;
 
		// kolejny poziom czasu aby zapisać wyniki (inaczej: ilosc minionego czasu z perspektywy animacji)
		double gifTime;

		//ilosc minionych iteracji
		int iters;

		// ilosc zapisanych wynikow (klatek animacji)	
		int frames;
		
		// dokladnosc w relaksacji cisnienia
		double EpsP;

		// maksymalna roznica cisnien P, zmienna uzywana w funkcji iteration w fragmencie ciśnienia
		double DiffPMax; 

		// tablica zapisujaca czy sasiedzi komorki sa pusci (uzywane w warunkach ciecz - atmosfera oraz funkcja emptyNeighbors)
		bool neighbors[4]; // N E S W 

		// suma pustych sasiadow komorki (uzywane w warunkach ciecz - atmosfera oraz funkcja emptyNeighbors)
		unsigned int neighSum;

		// tablice 

		//pozioma skladowa predkosci
		Matrix<double> U1;

		//pionowa skladowa predkosci
		Matrix<double> V1;

		//cisnienia 
		Matrix<double> P1;

		//niescisliwosc/dywergencja
		Matrix<double> D1;

		// flagi komorek F
		Matrix<Cell> F1;

		// żródła Z
		Matrix<double> Z1;

		// czastki znaczone C
		std::vector<Particle> C1;


		// dokladnosc zera uzywana przy przesuwaniu czastek znaczonych
		double zeroEps = 1e-08;
	
		// zmienna zmieniajaca metode liczenia dt na ekstrapolacje Richardsona
		bool dtER;

		// predkosc na wieczku (domyślnie 0.0)
		double lidSpeed;

		bool sU;
		bool sV;
		bool sP;
		bool sD;


    public:
		// główny konstruktor
		System(double w, double h, int x, int y, double gx, double gy, double v);

		// destruktor 
		~System();


		// dodanie brzegu - ramki dookoła
		void addBorder(bool slip); 

		//dodanie czastek znaczonych n - ilosc n x n 
		void addParticles(int n);

		// testowa scianka wody na start (argumenty są liczone w przedziałach siatki, nie wysokości/szerokości!)
		void fluidWall(int left, int right, int down, int up, int color=0); 

		// dodanie scianki - przeszkody w systemie - prostokat od left do right i down do up
		void solidWall(int left, int right, int down, int up, bool slip);

		//testowa kropla wody
		void fluidCircle(int x, int y, int r, int color=0);

		// dodanie wyjscia na dolnej sciance systemu na indeksie pos
		void addOut(int pos);

		// ustawienie prędkości na wieczku pojemnika
		void setLidSpeed(double s);

		// tworzenie pochylej cieczy - obszar lewy,dolny popd przekatna prostokatna wyznaczonego przez left, right, up, down
		// n x n czastek na komorke
		void fluidSlope(int left, int right, int down, int up, int n, int color=0);

		// stworzenie babelka powietrza w cieczy
		void bubble(int x, int y, int r, int color=0);


		// zapis do pliku flag/typów komórek
		void saveF();

		//zapis do pliku czastek
		void saveC() const;

		// generyczna funkcja zapisu tablicy Matrix do pliku o nazwie filename
		void saveMatrix(Matrix<double> m, std::string filename);

		// zapis do pliku prędkośći poziomych
		void saveU();

		// zapis do pliku prędkości pionowych
		void saveV();

		// zapis do pliku ciśnień
		void saveP();

		// zapis do pliku dywergencji
		void saveD();

		// zapis scianek/przeszkod do pliku
		void saveMap();


		// funkcja wypisujaca pomocnicze wartosci
        void print() const;

		// jedna iteracja o kroku czasowym dt
		void iteration(double dt, Matrix<double> &U, Matrix<double> &V, Matrix<double> &P,
			Matrix<double> &D, Matrix<Cell> &F, Matrix<double> &Z, std::vector<Particle> &C);

		// funkcja obsługująca warunki brzegowe dla ciśnień
		void boundaryConditionsP(int i, int j, Matrix<Cell> &F, Matrix<double> &P, Matrix<double> &U, Matrix<double> &V);  // WarunkiBrzegoweP

		// funkcja obsługująca pustych sąsiadów przy ropatrywaniu warunków ciecz - atmosfera
		void emptyNeighbors(int i, int j, Matrix<Cell> &F);

		// funkcja zwraca czy wszyscy sąsiedzi komórki o indeksach i j są pełni
		bool fullNeighbors(Matrix<int> &Is_any_c, int i, int j);

		// glowna funkcja przechodzaca przez krok w programie - iteracja + nowe dt
		void step();

		// krok w programie liczacy nowe dt z kryterium CFL
		void stepCFL();

		// krok w programie liczacy nowe dt z ekstrapolacji Richardsona
		void stepRE();

		// liczenie nowego kroku czasowego bazujace na estrapolacji Richardsona
		bool newDtRE(Matrix<double> U2, Matrix<double> V2);

		// liczenie nowego kroku czasowego bazujaca na kryterium CFL
		void newDtCFL();

		// funkcja pomocniczna, zwraca czy komorka po prawej gornej stronie obecnej jest brzegowa
		bool URON(Matrix<Cell> F, int i, int j); 

		// funkcja pomocnicza, zwraca czy ktorykolwiek sasiad komorki sasiaduje bezposrednio z brzegowa
		bool OB(Matrix<Cell> F, int i, int j); 

		// funkcja pomocnicza, zwraca czy komorka jest naroznikiem
		bool COR(Matrix<Cell> F, int i, int j); 

		// ustawienie liczenia dt na ekstrapolacje Richardsona
		void setER();

		// ustawienie liczenia dt na CFL (opcja domyslna)
		void setCFL();

		// zwraca miniony czas w symulacji
		double getTimePassed();

		void saveUVPD(bool u, bool v, bool p, bool d);

};

