#include "pch.h"
#include "Integration_Scheme_Adaptive.h"

namespace Com_Methods
{
	//�����������: �� ���� ������� ��� ������������ �������
	Integration_Scheme_Adaptive::Integration_Scheme_Adaptive(Integration_Scheme_Type Type)
	{
		//���������� �������� ����� � ����� ��������������
		switch (Type)
		{
			//����� ������ �����-1
		case Gauss1:
		{
			Weight = { 2 };
			Points = { Point(0, 0, 0) };
			break;
		}
		case Gauss2:
		{
			Weight = { 1,1 };
			Points = { Point(-1.0 / sqrt(3.0),0,0), Point(1.0 / sqrt(3.0),0,0) };
			break;
		}
		case Gauss3:
		{
			Weight = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
			Points = { Point(-sqrt(3.0 / 5.0),0,0),
					   Point(0,0,0),
					   Point(sqrt(3.0 / 5.0),0,0) };
			break;
		}
		case Simpson:
		{
			Weight = { 1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0 };
			Points = { Point(-1.0, 0, 0),
					   Point(0, 0, 0),
					   Point(1.0, 0, 0), };
			break;
		}
		}
	}

	//����� ��� ���������� ������������ ��������� �� ���������� �����: 
	//Begin � End - ������ � ����� ������� 
	//Num_Segments - ����� ���������
	//Func - ��������������� �������
	double Integration_Scheme_Adaptive::Calculate_Integral(
		const Point& Begin,
		const Point& End,
		int Number_Segments,
		double r,
		const std::function<double(const Point& P)>& Func) const
	{
		//��������� (������������ �����)
		double Result = 0.0;
		//��������� ����� ��������
		double X01,X02;
		//����� �����, ������������ ����� � ������ �������
		double power = 1, sum1 = 0, sum2 = 0;
		for (int i = 0; i < Number_Segments/2; i++)
		{
			sum1 += power;
			sum2 += 1.0 / power;
			power *= r;
		}
		double h1 = (End.x() - Begin.x())/ 2.0 / sum1, 
			h2 = (End.x() - Begin.x())/ 2.0 / sum2;
		//����� �� ���� ��������� ���������
		power = 1;
		X01 = Begin.x();
		X02 = (End.x() - Begin.x()) / 2.0;
		for (int Integ_Point = 0; Integ_Point < Points.size(); Integ_Point++)
		{
			//������� � ������-�������� [-1, 1]
			auto P1 = Point(X01 + (1 + Points[Integ_Point].x()) * h1 * power / 2.0, 0, 0);
			auto P2 = Point(X02 + (1 + Points[Integ_Point].x()) * h2 * (1.0 / power) / 2.0, 0, 0);
			Result += Weight[Integ_Point] * Func(P1)* (h1*power / 2.0);
			Result += Weight[Integ_Point] * Func(P2)* (h2*(1.0/power) / 2.0);
		}
		for (int i = 1; i < Number_Segments/2; i++)
		{
			//��������� ����� ��������
			X01 += h1*power;
			X02 += h2 * (1 / power);
			//����� �� ����� ��������������
			for (int Integ_Point = 0; Integ_Point < Points.size(); Integ_Point++)
			{
				//������� � ������-�������� [-1, 1]
				auto P1 = Point(X01 + (1 + Points[Integ_Point].x()) * h1*power / 2.0, 0, 0);
				auto P2 = Point(X02 + (1 + Points[Integ_Point].x()) * h2 * (1.0/power) / 2.0, 0, 0);
				Result += Weight[Integ_Point] * Func(P1) * (h1 * power / 2.0);
				Result += Weight[Integ_Point] * Func(P2) * (h2 * (1.0 / power) / 2.0);
			}
			power *= r;
		}
		//��������� ��������� � ������ �������� �� ������� [-1, 1]
		return Result;
	}
}
