
#include <iostream>
#include "fluid.h"
#include "scene.h"
#include "utilities.h"
#include "string.h"
/*
	Глобальные переменные
*/
using enum debug_mode;

unsigned int num_cells_x = 100; // размер решётки по горизонтали
unsigned int num_cells_y = 100; // размер решётки по вертикали

/*
	Конец глобальных переменных
*/

/* 
    input:    	Аргументы командной строки

    output:   	нет.

    remarks:	Обрабатывает аргументы командной строки для заполнения глобальных переменных
*/
void process_arguments(int argc, char const *argv[])
{
	for (int i {}; i < argc; i++)
	{
		if (strcmp(argv[i], "-x") == 0)
		{
			num_cells_x = std::stoi(argv[i+1]);
		}
		if (strcmp(argv[i], "-y") == 0)
		{
			num_cells_y = std::stoi(argv[i+1]);
		}
	}
}

void single_var_printout(std::string name, double var)
{
	std::cout << name << ": " << var << std::endl;
}

void config_printout()
{
	std::cout << "Программа симуляции будет запущена со следующими параметрами:" << std::endl;
	single_var_printout("Размер сетки по оси Х", num_cells_x);
	single_var_printout("Размер сетки по оси Х", num_cells_y);
	
}

int main(int argc, char const *argv[])
{	
	process_arguments(argc, argv); // обработка аргументов
	config_printout();
	Scene test_scene(num_cells_x,num_cells_y);
	test_scene.add_boundary_inlet(5);
	for (int i = 0; i < 60; i++)
	{
		test_scene.debug_print();
		test_scene.simulate();
	}
	
	return 0;
}