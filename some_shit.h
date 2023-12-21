#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>


using namespace std;
namespace origin{
	constexpr double c2 = 0.05;
	constexpr double A = 3;
	constexpr double B = 3;
	constexpr double C = -3;

	//if you want to change this -> be careful
	//with runge_kutta_2p_2e!!!
	constexpr double x0 = 0;
	constexpr double x1 = 5;
	const vector<double> y0 = { 1, 1, A, 1 };
};

void print(const vector<double>& v) {
	for (auto i : v) cout << i << " ";
}
vector<double> root_sol(double x) {
	return { exp(sin(x * x))
		, exp(origin::B*sin(x*x))
		, origin::C * sin(x * x) + origin::A
		, cos(x * x) };
}

vector<double> f(double x, const vector<double>& y) {
	return { 2 * x * pow(y[1], 1.0 / origin::B) * y[3]
	, 2 * origin::B * x * exp(origin::B / origin::C * (y[2] - origin::A)) * y[3]
	, 2 * origin::C * x * y[3]
	, -2 * x * log(y[0]) };
}
double part_f1(double x, const vector<double>& y) {
	return 2 * x * pow(y[1], 1.0 / origin::B) * y[3];
}
double part_f2(double x, const vector<double>& y) {
	return 2 * origin::B * x * exp(origin::B / origin::C * (y[2] - origin::A)) * y[3];
}
double part_f3(double x, const vector<double>& y) {
	return 2 * origin::C * x * y[3];
}
double part_f4(double x, const vector<double>& y) {
	return -2 * x * log(y[0]);
}

//to do
void output_to_file() {
	return;
}
void runge_kutta_2(double h) {
	ofstream fout("data.txt");
	ofstream fout_or("data_or.txt");
	fout << "{";
	fout_or << "{";
	vector<double> y_n = origin::y0;
	double a21 = origin::c2;
	double b2 = 1.0 / (2.0 * origin::c2);
	double b1 = 1 - 1.0 / (2.0 * origin::c2);

	cout << "params:\n";
	cout << a21 << " " << b2 << " " << b1 << "\n";
	for (double x = 0; x <= origin::x1; x += h) {
		fout << "'" << x << "':[";
		fout_or << "'" << x << "':[";
		//here i need to compute y_(n+1) based on y_n
		
		vector<double> K1 = f(x, y_n);
		vector<double> helper = y_n;
		helper[0] += h * a21 * K1[0];
		helper[1] += h * a21 * K1[1];
		helper[2] += h * a21 * K1[2];
		helper[3] += h * a21 * K1[3];

		vector<double> K2 = f(x, helper);
		vector or_val = root_sol(x);
		for (int el = 0; el < 4; el++) {
			if (el != 3) {
				fout << y_n[el] << ",";
				fout_or << or_val[el] << ",";
			}
			else {
				fout << y_n[el] << "],";
				fout_or << or_val[el] << "],";	
			}
		}
		y_n[0] = y_n[0] + h * b1 * K1[0] + h * b2 * K2[0];
		y_n[1] = y_n[1] + h * b1 * K1[1] + h * b2 * K2[1];
		y_n[2] = y_n[2] + h * b1 * K1[2] + h * b2 * K2[2];
		y_n[3] = y_n[3] + h * b1 * K1[3] + h * b2 * K2[3];
	}
	fout << "}";
	fout_or << "}";
	fout.close();
	fout_or.close();
}

void runge_kutta_4(double h) {
	ofstream fout("data.txt");
	ofstream fout_or("data_or.txt");
	fout << "{";
	fout_or << "{";
	vector<double> y_n = origin::y0;
	double a21 = origin::c2;
	double b2 = 1.0 / (2.0 * origin::c2);
	double b1 = 1 - 1.0 / (2.0 * origin::c2);

	cout << "params:\n";
	cout << a21 << " " << b2 << " " << b1 << "\n";
	for (double x = 0; x <= origin::x1; x += h) {
		fout << "'" << x << "':[";
		fout_or << "'" << x << "':[";
		//here i need to compute y_(n+1) based on y_n

		vector<double> K1 = f(x, y_n);
		vector<double> helper = y_n;
		helper[0] += h * a21 * K1[0];
		helper[1] += h * a21 * K1[1];
		helper[2] += h * a21 * K1[2];
		helper[3] += h * a21 * K1[3];

		vector<double> K2 = f(x, helper);
		vector or_val = root_sol(x);
		for (int el = 0; el < 4; el++) {
			if (el != 3) {
				fout << y_n[el] << ",";
				fout_or << or_val[el] << ",";
			}
			else {
				fout << y_n[el] << "],";
				fout_or << or_val[el] << "],";
			}
		}
		y_n[0] = y_n[0] + h * b1 * K1[0] + h * b2 * K2[0];
		y_n[1] = y_n[1] + h * b1 * K1[1] + h * b2 * K2[1];
		y_n[2] = y_n[2] + h * b1 * K1[2] + h * b2 * K2[2];
		y_n[3] = y_n[3] + h * b1 * K1[3] + h * b2 * K2[3];
	}
	fout << "}";
	fout_or << "}";
	fout.close();
	fout_or.close();
}