#include "Cell.h"
#include <iostream>

Cell::Cell(){
	type = emptyc;
}

Cell::Cell(Type t):type(t){}

void Cell::set(Type t){
	type = t;
}

Type Cell::get() const {
	return type;
}

void Cell::setColor(int c)
{
	color = c;
};

int Cell::getColor() const 
{
	return color;
};



bool Cell::operator == (const Cell &k)
{
	if (type == k.get()){
		return true;
	}
	return false;
}

bool Cell::operator != (const Cell &k)
{
	if (type != k.get()){
		return true;
	}
	return false;
}

bool Cell::boundary() const{
	return (type == border || type == border_noslip  || type == out);
};

bool Cell::slip() const{
	return (type == border);
};

bool Cell::noSlip() const{
	return (type == border_noslip);
};

bool Cell::fluid() const{
	return (type == full || type == surface);
};

