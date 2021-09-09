#include "pch.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <functional>
#include "Point.h"
#include "Integration_Scheme_Interval.h"
#include <string>

int main()
{
	//путь записи
	const std::string PATH = "results/";
	//имя функции
	const std::string f_name = "ex_";
	//подынтегральная функция f(x) = e^x
	std::function<double(const Com_Methods::Point &P)> f = 
	[](const Com_Methods::Point &P) { return std::exp(P.x()); };
	//первообразная F(x) = e^x
	std::function<double(const Com_Methods::Point &P)> F =
	[](const Com_Methods::Point &P) { return std::exp(P.x()); };


	//квадратурная формула Гаусс-1
	Com_Methods::Integration_Scheme_Interval Gauss1(Com_Methods::Integration_Scheme::Gauss1);
	//квадратурная формула Гаусс-2
	Com_Methods::Integration_Scheme_Interval Gauss2(Com_Methods::Integration_Scheme::Gauss2);
	//квадратурная формула Гаусс-3
	Com_Methods::Integration_Scheme_Interval Gauss3(Com_Methods::Integration_Scheme::Gauss3);
	//квадратурная формула Симпсон
	Com_Methods::Integration_Scheme_Interval Simpson(Com_Methods::Integration_Scheme::Simpson);

	//начало и конец отрезка интегрирования
	auto Begin = Com_Methods::Point(0, 0, 0);
	auto End   = Com_Methods::Point(1, 0, 0);

	//число сегментов
	int Num_Segments = 10;
	
	//точное значение интеграла (ф. Ньютона-Лейбница)
	double I_True = F(End) - F(Begin);

	//численное значение интеграла
	double I[3];

	//Нахождение интегралом методом Гаусс-1
	Num_Segments = 10;
	std::ofstream Writer(PATH+f_name+"Gauss1.txt");
	Writer.setf(std::ios::scientific);
	Writer << "h\tI\t|I - I_True|"<< std::endl;
	int i = 0;
	I[i] = Gauss1.Calculate_Integral(Begin, End, Num_Segments, f);

	Writer << (End.x() - Begin.x()) / Num_Segments << "\t";
	Writer << I[i] << "\t";
	Writer << fabs(I[i] - I_True) << std::endl;
	Num_Segments *= 2;
	i++;
	for (; i < 3; i++)
	{
		I[i] = Gauss1.Calculate_Integral(Begin, End, Num_Segments, f);

		Writer << (End.x() - Begin.x()) / Num_Segments << "\t";
		Writer << I[i] << "\t";
		Writer << fabs(I[i] - I_True) << "\t";
		Writer << fabs(I[i-1] - I_True) / fabs(I[i] - I_True) << std::endl;
		Num_Segments *= 2;
	}
	Writer.close();
	
	//Нахождение интегралом методом Гаусс-2
	Num_Segments = 10;
	Writer.open(PATH + f_name + "Gauss2.txt");
	Writer.setf(std::ios::scientific);
	Writer << "h\tI\t|I - I_True|" << std::endl;
	i = 0;
	I[i] = Gauss2.Calculate_Integral(Begin, End, Num_Segments, f);

	Writer << (End.x() - Begin.x()) / Num_Segments << "\t";
	Writer << I[i] << "\t";
	Writer << fabs(I[i] - I_True) << std::endl;
	Num_Segments *= 2;
	i++;
	for (; i < 3; i++)
	{
		I[i] = Gauss2.Calculate_Integral(Begin, End, Num_Segments, f);

		Writer << (End.x() - Begin.x()) / Num_Segments << "\t";
		Writer << I[i] << "\t";
		Writer << fabs(I[i] - I_True) << "\t";
		Writer << fabs(I[i - 1] - I_True) / fabs(I[i] - I_True) << std::endl;
		Num_Segments *= 2;
	}
	Writer.close();

	//Нахождение интегралом методом Гаусс-3
	//первый вариант для Гаусс-3
	Num_Segments = 10;
	Writer.open(PATH + f_name + "Gauss3.txt");

	////вариант с более мелким шагом для Гаусс-3
	//Num_Segments = 40;
	//Writer.open(PATH + f_name + "Gauss3_1.txt");

	Writer.setf(std::ios::scientific);
	Writer << "h\tI\t|I - I_True|" << std::endl;
	i = 0;
	I[i] = Gauss3.Calculate_Integral(Begin, End, Num_Segments, f);

	Writer << (End.x() - Begin.x()) / Num_Segments << "\t";
	Writer << I[i] << "\t";
	Writer << fabs(I[i] - I_True) << std::endl;
	Num_Segments *= 2;
	i++;
	for (; i < 3; i++)
	{
		I[i] = Gauss3.Calculate_Integral(Begin, End, Num_Segments, f);

		Writer << (End.x() - Begin.x()) / Num_Segments << "\t";
		Writer << I[i] << "\t";
		Writer << fabs(I[i] - I_True) << "\t";
		Writer << fabs(I[i - 1] - I_True) / fabs(I[i] - I_True) << std::endl;
		Num_Segments *= 2;
	}
	Writer.close();

	//Нахождение интегралом методом Симпсона
	Num_Segments = 10;
	Writer.open(PATH + f_name + "Simpson.txt");
	Writer.setf(std::ios::scientific);
	Writer << "h\tI\t|I - I_True|" << std::endl;
	i = 0;
	I[i] = Simpson.Calculate_Integral(Begin, End, Num_Segments, f);

	Writer << (End.x() - Begin.x()) / Num_Segments << "\t";
	Writer << I[i] << "\t";
	Writer << fabs(I[i] - I_True) << std::endl;
	Num_Segments *= 2;
	i++;
	for (; i < 3; i++)
	{
		I[i] = Simpson.Calculate_Integral(Begin, End, Num_Segments, f);

		Writer << (End.x() - Begin.x()) / Num_Segments << "\t";
		Writer << I[i] << "\t";
		Writer << fabs(I[i] - I_True) << "\t";
		Writer << fabs(I[i - 1] - I_True) / fabs(I[i] - I_True) << std::endl;
		Num_Segments *= 2;
	}
	Writer.close();


	//Исследование порядка точности
	//полином вида (m+1)*x^m
	//подынтегральная функция f(x) = e^x
	int m = 1;
	std::function<double(const Com_Methods::Point& P)> f_polinom =
		[&](const Com_Methods::Point& P) { return (m+1)*std::pow(P.x(),m); };
	//первообразная F(x) = e^x
	std::function<double(const Com_Methods::Point& P)> F_polinom =
		[&](const Com_Methods::Point& P) { return std::pow(P.x(), m+1); };
	
	//функция исследования и вывода результатов в файл
	std::function<
		void(Com_Methods::Integration_Scheme_Interval &Integral, const Com_Methods::Point& a, const Com_Methods::Point& b, const int M, int NumSegments, const std::string &path)> Output3=
		[&](Com_Methods::Integration_Scheme_Interval& Integral, const Com_Methods::Point& a, const Com_Methods::Point& b, const int M, int NumSegments, const std::string& path) {

		Writer.open(path);
		Writer << "Степень полинома\t";
		Writer << "f(x)\t";
		Writer << "Шаг h\t";
		Writer << "I*-Ih\t";
		Writer << "(I*-Ih)/(I*-Ih/2)\t";
		Writer << "(Ih/2-Ih)/(2^k-1)\t";
		Writer << "I^R\t";
		Writer << "I*-I^R" << std::endl;
		for (m = M - 1; m <= M + 1; m++)
		{
			//Точное значение интеграла
			I_True = F_polinom(b) - F_polinom(a);
			//Шаг h
			Num_Segments = NumSegments;
			for (i=0; i < 3; i++)
			{
				I[i] = Integral.Calculate_Integral(a, b, Num_Segments, f_polinom);
				Num_Segments *= 2;
			}
			Num_Segments = NumSegments;
			for (i = 0; i < 2; i++)
			{
				double IR;
				Writer << m << "\t";//степень полинома
				Writer << std::to_string(m+1) + "x^" + std::to_string(m) + "\t";//f(x)
				Writer << std::to_string((b.x() - a.x()) / Num_Segments) + "\t";//Шаг h
				Writer << I_True - I[i] << "\t";//I*-Ih
				Writer << (I_True - I[i])/(I_True-I[i+1]) << "\t";//(I*-Ih)/(I*-Ih/2)
				IR = (I[i + 1] - I[i]) / (std::pow(2, M+1) - 1);
				Writer << IR << "\t";//(Ih/2-Ih)/(2^k-1)
				IR += I[i + 1];
				Writer << IR << "\t";//IR
				Writer << I_True - IR << "\t";//I_true-TR
				Writer << std::endl;
				Num_Segments *= 2;
			}
			Writer << m << "\t";//степень полинома
			Writer << std::to_string(m+1) + "x^" + std::to_string(m) + "\t";//f(x)
			Writer << std::to_string((b.x() - a.x()) / Num_Segments) + "\t";//Шаг h
			Writer << I_True - I[i] << "\t";//I*-Ih
			Writer << "-"<< "\t";//(I*-Ih)/(I*-Ih/2)
			Writer << "-" << "\t";//(Ih/2-Ih)/(2^k-1)
			Writer << "-" << "\t";//IR
			Writer << "-" << "\t";//I_true-TR
			Writer << std::endl;
		}
		Writer.close();

	};
	int NumSegments = 10;
	Com_Methods::Point a(0,0,0), b(1,0,0);
	Output3(Gauss1, a,b, 1, NumSegments, PATH + "Gauss1.txt");
	Output3(Gauss2, a,b, 3, NumSegments, PATH + "Gauss2.txt");
	Output3(Gauss3, a,b, 5, NumSegments, PATH + "Gauss3.txt");
	Output3(Simpson, a,b, 3, NumSegments, PATH + "Simpson.txt");

	//Иследования неполиноминальной функции xsin(10000x)
	std::function<double(const Com_Methods::Point& P)> f_xsin =
		[](const Com_Methods::Point& P) {return P.x() * std::sin(10000 * P.x()); };
	std::function<double(const Com_Methods::Point& P)> F_xsin =
		[](const Com_Methods::Point& P) {return std::sin(10000 * P.x())/100000000 - P.x()*std::cos(10000*P.x())/10000; };
	//функция вывода результатов
	std::function<
		void(Com_Methods::Integration_Scheme_Interval& Integral, const Com_Methods::Point& a, const Com_Methods::Point& b, const int M, int NumSegments, const std::string& path)> Output4 =
		[&](Com_Methods::Integration_Scheme_Interval &Integral, const Com_Methods::Point& a, const Com_Methods::Point& b, const int M, int NumSegments, const std::string& path) {
		Writer.open(path);
		Writer << "Шаг h\t";
		Writer << "I*-Ih\t";
		Writer << "(I*-Ih)/(I*-Ih/2)\t";
		Writer << "(Ih/2-Ih)/(2^k-1)\t";
		Writer << "I^R\t";
		Writer << "I*-I^R" << std::endl;
		//Точное значение интеграла
		I_True = F_xsin(b) - F_xsin(a);
		//Шаг h
		Num_Segments = NumSegments;
		for (i = 0; i < 3; i++)
		{
			I[i] = Integral.Calculate_Integral(a, b, Num_Segments, f_xsin);
			Num_Segments *= 2;
		}
		Num_Segments = NumSegments;
		for (i = 0; i < 2; i++)
		{
			double IR;
			Writer << std::to_string((b.x() - a.x()) / Num_Segments) + "\t";//Шаг h
			Writer << I_True - I[i] << "\t";//I*-Ih
			Writer << (I_True - I[i]) / (I_True - I[i + 1]) << "\t";//(I*-Ih)/(I*-Ih/2)
			IR = (I[i + 1] - I[i]) / (std::pow(2, M + 1) - 1);
			Writer << IR << "\t";//(Ih/2-Ih)/(2^k-1)
			IR += I[i + 1];
			Writer << IR << "\t";//IR
			Writer << I_True - IR << "\t";//I_true-TR
			Writer << std::endl;
			Num_Segments *= 2;
		}
		Writer << std::to_string((b.x() - a.x()) / Num_Segments) + "\t";//Шаг h
		Writer << I_True - I[i] << "\t";//I*-Ih
		Writer << "-" << "\t";//(I*-Ih)/(I*-Ih/2)
		Writer << "-" << "\t";//(Ih/2-Ih)/(2^k-1)
		Writer << "-" << "\t";//IR
		Writer << "-" << "\t";//I_true-TR
		Writer << std::endl;
		Writer.close();
	};
	NumSegments = 10000;
	Output4(Gauss1, a, b, 1, NumSegments, PATH + "sinx_nested_Gauss1.txt");
	Output4(Gauss2, a, b, 3, NumSegments, PATH + "sinx_nested_Gauss2.txt");
	Output4(Gauss3, a, b, 5, NumSegments, PATH + "sinx_nested_Gauss3.txt");
	Output4(Simpson, a, b, 3, NumSegments, PATH + "sinx_nested_Simpson.txt");
	Writer.open(PATH + "xsin_adaptive.txt");
	I_True = F_xsin(b) - F_xsin(a);
	Writer << "I*\tI\t|I*-I|\t|I*-I|/I*"<<std::endl;
	double r = 1.0;
	I[0] = Simpson.Calculate_Integral_Adaptive(a, b, NumSegments, r, f_xsin);
	Writer << I_True << "\t " << I[0] << "\t" << r << "\t" << std::abs(I_True-I[0]) << "\t" << (std::abs((I_True - I[0])/I_True))<<std::endl;
	//коэфициент разрядки
	r = 1.4;
	I[0] = Simpson.Calculate_Integral_Adaptive(a, b, NumSegments, r, f_xsin);
	Writer << I_True << "\t " << I[0] << "\t" << r << "\t" << std::abs(I_True - I[0]) << "\t" << (std::abs((I_True - I[0]) / I_True)) << std::endl;
	r = 1.5;
	I[0] = Simpson.Calculate_Integral_Adaptive(a, b, NumSegments, r, f_xsin);
	Writer << I_True << "\t " << I[0] << "\t" << r << "\t" << std::abs(I_True - I[0]) << "\t" << (std::abs((I_True - I[0]) / I_True)) << std::endl;
	r = 1.6;
	I[0] = Simpson.Calculate_Integral_Adaptive(a, b, NumSegments, r, f_xsin);
	Writer << I_True << "\t " << I[0] << "\t" << r << "\t" << std::abs(I_True - I[0]) << "\t" << (std::abs((I_True - I[0]) / I_True)) << std::endl;
	Writer.close();
	return 0;
}