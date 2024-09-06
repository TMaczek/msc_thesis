#pragma once
#include <iostream>
#include "Type.h"

/**
 * @brief 
 * Klasa reprezentująca komórkę na siatce Eulerowskiej. 
 */
class Cell{
	private:
		// typ komórki
		Type type;
		// kolor komórki - używany przy dodawaniu obszarów cieczy
		int color;

	public:		
		Cell();

		Cell(Type t);

		void set(Type t);

		Type get() const;

		void setColor(int c);

		int getColor() const;

		// operator przyrównujący do typu
		bool operator == (const Cell &k);

		// operator przyrównujący do typu
		bool operator != (const Cell &k);

		// funkcja sprawdzająca czy komórka jest jedną z brzegowych
		bool boundary() const;

		// funkcja sprawdzająca czy komórka jest z poślizgiem
		bool slip() const; 

		// funkcja sprawdzająca czy komórka jest bez poślizgu
		bool noSlip() const;
		
		// funkcja sprawdzająca czy komórka zawiera ciecz
		bool fluid() const; 

};