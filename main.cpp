#include <iostream>

#include "add.h"

int main()
{
	const int arraySize = 5;
	const int a[arraySize] = { 1, 2, 3, 4, 5 };
	const int b[arraySize] = { 10, 20, 30, 40, 50 };
	int c[arraySize] = { 0 };

	add(c, a, b, arraySize);

	for (int i = 0; i < arraySize; i++)
	{
		std::cout << c[i] << std::endl;
	}

	return 0;
}