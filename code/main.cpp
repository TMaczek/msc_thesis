#include <iostream>
#include "System.h"
#include <cmath>


using namespace std;

int main(int argc, char *argv[]){
	// przykladowe konfiguracje
    
	// 1) Opuszczona tama
	/*

	System cont1(200, 200, 100, 100, 0, -0.1, 1.002); 
	cont1.addBorder(false); 	//true = poslizg
	cont1.fluidWall(1, 30, 1, 50);
	cont1.addParticles(2); // 4 x 4 = 16 czastek na komorke cieczy
 	cont1.saveMap(); // zapis mapki do pliku
	for(int i=0; i<10000; ++i) // mozna wg ilosci iteracji lub wykorzystac getTimePassed
	{
		cont1.step();
		if(i%500 == 0)
			cout << "Iteracja: " << i << endl;
	}	

	*/

	// 2) Spadająca kropla
	// 2_1 - głęboka (depth 60), 2_2 - płytka (depth 15), r=8
	 /*
	System cont2(200, 200, 100, 100, 0, -0.1, 1.002);
	cont2.addBorder(true);
	cont2.fluidWall(1, 99, 1, 30, 1); // usunac dla jednej wersji
	cont2.fluidCircle(50, 70, 7, 2);
	cont2.addParticles(2);
	cont2.saveMap();
	for(int i=0; i<2000; ++i)
	{
		cont2.step();
		if(i%500 == 0)
			cout << "Iteracja: " << i << endl;

	}
	 */

	// 3) Wyrównywanie się poziomów
	/*
	System cont3(200, 200, 100, 100, 0, -0.1, 1.002);
	cont3.addBorder(true);
	cont3.solidWall(23, 27, 20, 99, true);
	cont3.solidWall(48, 52, 20, 99, true);
	cont3.solidWall(73, 77, 20, 99, true);
	
	cont3.fluidWall(1, 99, 1, 20, 0);
	cont3.fluidWall(1, 23, 20, 70, 1);
	cont3.fluidWall(27, 48, 20, 50, 2);
	cont3.fluidWall(52, 73, 20, 30, 3);
	cont3.fluidWall(77, 99, 20, 90, 4);

	cont3.addParticles(1);
	cont3.saveMap();
	cont3.saveUVPD(false, false, true, false); // zapis cisnien do pliku
	for(int i=0; i<8000; ++i)
	{
		cont3.step();
		if(i%500 == 0)
			cout << "Iteracja: " << i << endl;
	}
	*/

	// 4) Rozbijanie się fal
	/*
	double g_abs = 0.1;
	double gy = -0.5 * std::pow(3*g_abs, 0.5) ;
	double gx = -0.5 * std::pow(g_abs, 0.5);

	System cont4(500, 300, 100, 60, gx, gy, 1.002);
	cont4.addBorder(false); 
	cont4.fluidWall(1, 20, 1, 45);
	cont4.fluidSlope(20, 62, 1, 25, 3);	
	cont4.addParticles(3);
	cont4.saveMap();
	for(int i=0; i<6000; ++i)
	{
		cont4.step();
	}
	cont4.saveF();
	*/

	// 5) Unoszący się bąbelek
	/*
	System cont5(50, 200, 25, 100, 0.0, -0.1, 1.002);
	cont5.addBorder(true);
	cont5.bubble(1, 11, 8, 1);
	cont5.addParticles(2);
	cont5.saveMap();
	for(int i=0; i<2000; ++i)
	{
		cont5.step();
		if(i%200==0)
			cout << i << endl;

	}
	*/


	// 6) Pojemnik z prędkością na wieczku
	// zapisywac predkosci i cisnienia?
	/*
	System cont6(100, 100, 50, 50, 0, -0.1, 1.002); // siatka wczesniej byla 40 x 40
	cont6.addBorder(false);
	cont6.fluidWall(1, 49, 1, 49);
	cont6.setLidSpeed(4.0);
	cont6.addParticles(1);
	cont6.saveMap();
	cont6.saveUVPD(true, true, false, false);
	for(int i=0; i<7000; ++i)
	{
		cont6.step();
		if(i%500 == 0)
			cout << "Iteracja: " << i << endl;
	}
	*/

	// 7) Przeszkoda i stopień 
	/*
	System cont7(300, 200, 150, 100, 0, -0.1, 2.12);  // lepkość 2.12 lub 1.002
	cont7.solidWall(75, 80, 1, 6, true);
	cont7.fluidWall(1, 30, 1, 50);
	cont7.addParticles(2); // 4 x 4 = 16 czastek na komorke cieczy
 	cont7.saveMap(); // zapis mapki do pliku
	for(int i=0; i<5000; ++i) 
	{
		cont7.step();
	}
	*/

	// rozbijanie fali
	/*
	System cont7b(300, 200, 150, 100, 0, -0.1, 1.002); 
	cont7b.addBorder(true); 	//true = poslizg
	cont7b.solidWall(30, 35, 20, 99, true);
	cont7b.fluidWall(1, 30, 1, 95, 0);
	cont7b.fluidWall(30, 149, 1, 11, 4); // 18 lub 4 lub 11/12
	cont7b.addParticles(2);
	cont7b.saveMap();

	for(int i=0; i<3500; ++i)
	{
		cont7b.step();
	}
	*/

	// 8) Output 
	/*
	System cont8(200, 100, 100, 50, 0, -0.1, 1.002); // siatka wczesniej byla 40 x 40
	cont8.addBorder(true);
	cont8.addOut(48);
	cont8.addOut(49);
	cont8.addOut(50);
	cont8.addOut(51);
	cont8.addOut(52);
	cont8.addOut(53);
	cont8.fluidWall(1, 99, 1, 20, 0); 
	cont8.addParticles(1); // za duza liczba czastek rowniez psuje symulacje
	cont8.saveMap();
	cont8.saveUVPD(false, false, true, false);
	for(int i=0; i<4000; ++i)
	{
		cont8.step();
		if(i%200==0)
			cout << i << endl;
	}
	*/

    return 0;
}
