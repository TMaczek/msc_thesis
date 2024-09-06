#pragma once
#include <string>

/**
 * @brief 
 * Klasa reprezentująca cząstkę znaczoną w programie.
 * W razie rozbudowy programu o wejście warto dodać tu flagę oznaczającą czy cząstka dalej jest w komórce wejściowej.
 * 
 */
class Particle{
    private:
        double x;		// Współrzędna x
		double y;		// Współrzędna y
		double uk; 		// Prędkość w kierunku osi x
		double vk; 		// Prędkość w kierunku osi y
		int color;		// Kolor cząstki

    public:
		Particle();

		Particle(double x, double y, double uk, double vk);

		Particle(double x, double y, double uk, double vk, int color);

		double getX() const;
		double getY() const;
		double getUk() const;
		double getVk() const;
		
		void setX(double xn);
		void setY(double yn);
		void setUk(double ukn);
		void setVk(double vkn);

		// wypisywanie danych cząstki, używane przy zapisie do pliku
		std::string printToFile() const;
};
