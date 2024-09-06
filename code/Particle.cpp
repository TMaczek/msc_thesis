#include "Particle.h"
#include <iostream>
#include <string>


Particle::Particle():x(0.0), y(0.0), uk(0.0), vk(0.0), color(1){};

Particle::Particle(double x, double y, double uk, double vk):
	x(x), y(y), uk(uk), vk(vk), color(1){};


Particle::Particle(double x, double y, double uk, double vk, int color):
	x(x), y(y), uk(uk), vk(vk), color(color){};


std::string Particle::printToFile() const{ 
	return std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(color);
}

double Particle::getX() const 
{
	return x;
};

double Particle::getY() const 
{
	return y;
};
		
double Particle::getUk() const 
{
	return uk;
};

double Particle::getVk() const 
{
	return vk;
};

void Particle::setX(double xn){
	x = xn;
}

void Particle::setY(double yn){
	y = yn;
}

void Particle::setUk(double ukn){
	uk = ukn;
}

void Particle::setVk(double vkn){
	vk = vkn;
}