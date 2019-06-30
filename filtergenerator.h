#pragma once

double GetXn(double x, double width);
double GetYn(double y, double height);
double Factorial(double x);

Mat GenerateHF0_f(int width, int height, double x1, double x2);			//to generate first HF(slightely different from other HF)
Mat GenerateHF_f(int width, int height, double x1, double x2);
Mat GenerateLF0_f(int width, int height, double x1, double x2);			//to generate first LF(slightely different from other LF)
Mat GenerateLF_f(int width, int height, double x1, double x2);
Mat GenerateOF_f(int width, int height, int orientInd, int orientMax);