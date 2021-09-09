#pragma once
#include "Integration_Scheme.h"
#include <functional>

namespace Com_Methods
{
	class Integration_Scheme_Adaptive : protected Integration_Scheme
	{
	public:
		//�����������: �� ���� ������� ��� ������������ �������
		Integration_Scheme_Adaptive(Integration_Scheme_Type Type);
		//����� ��� ���������� ������������ ��������� �� ���������� �����: 
		//Begin � End - ������ � ����� ������� 
		//Num_Segments - ����� ���������
		//Func - ��������������� �������
		//r - ����������� ��������
		double Calculate_Integral(const Point& Begin,
			const Point& End,
			int Number_Segments,
			double r,
			const std::function<double(const Point& P)>& Func) const;
	};
}

