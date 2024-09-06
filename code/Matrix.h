#pragma once
#include <iostream>
#include <vector>

// https://stackoverflow.com/questions/936687/how-do-i-declare-a-2d-array-in-c-using-new


template <typename T> class Matrix{
	private:
		std::vector<T> arr;

		int width;
		int height;
	
	public:
		Matrix(){
			width = 0;
			height = 0;

		}

		Matrix(int i, int j){
			init(i, j);
		}; 

		void init(int i, int j){
			if (arr.size()==0){
				arr = std::vector<T>(i * j);
				width = i;
				height = j;

			}
			else{ 
				std::cout << "Matrix: pamięć już zaalokowana!" << std::endl;
			}
		}

		int getWidth() const{ return width;}

		int getHeight() const {return height;}

		~Matrix(){
			// Dodać w razie potrzeby
		}; 


		void copy(Matrix<T> m){
			init(m.getWidth(), m.getHeight());
		}

		// funkcja indeksowania
		T at(int x, int y) const{ 
			return arr[x * width + y];
		};

		// ustawianie elementow
		void set(int x, int y, T val){
			arr[x + width * y] = val;
		}

		// zwiekszenie elementu
		void add(int x, int y, T val){
			arr[x + width * y] += val;
		}

		// operator indeksowania - wspolrzedne w nawiasach okraglych, np M(0, 1)
		T& operator()(int x, int y){
			return arr[x + height *y];
		}

		T* operator[](int i){
			return &(arr[i * height]); 
		}

};
