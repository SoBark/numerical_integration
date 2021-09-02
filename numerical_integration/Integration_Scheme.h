#pragma once
#include "Point.h"
#include <vector>

namespace Com_Methods
{
	//класс схемы интегрирования
	class Integration_Scheme
	{
	protected:
		//узлы интегрирования
		std::vector<Point> Points;
		//веса квадратурной формулы
		std::vector<double> Weight;
	public:
		//типы квадратурных формул
		enum Integration_Scheme_Type
		{
			Gauss1 = 1,
			Gauss2,
			Gauss3,
			Simpson
		};
	};
}

