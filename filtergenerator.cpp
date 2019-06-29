#include "stdafx.h"
#include "filtergenerator.h"

double GetXn(double x, double width)
{
	return x - width/2;
}

double GetYn(double y, double height)
{
	return -(y - height/2);
}

double Factorial(double x)
{
	double sum = 1.0;
	for (int i = 1; i <= x; ++i)
	{
		sum *= i;
	}
	return sum;
}

Mat GenerateHF0_f(int width, int height, double x1, double x2)
{
	Mat HF = Mat(width, height, CV_32FC1);

	for (double y = 0; y < height; y++)
	{
		for (double x = 0; x < width; x++)
		{
			double xn = GetXn(x, width);
			double yn = GetYn(y, height);

			double radius = sqrt(xn * xn + yn * yn);
			double ttmp2;

			if (radius <= x1) { ttmp2 = 0.0; }
			else if (radius >= x2) { ttmp2 = 1.0; }
			else // if ((radius > x1) && (radius < x2))
			{
				double ttmp = M_PI / 4.0 * (1.0 + (radius - x1) / (x2 - x1));
				ttmp2 = cos(M_PI / 2.0 * (log(ttmp * 2.0 / M_PI) / log(2.0)));
			}

			HF.at<float>(y, x) = ttmp2;
		}
	}

	Mat planes[] = { HF,HF };
	merge(planes, 2, HF);
	return HF;
}

Mat GenerateHF_f(int width, int height, double x1, double x2)
{
	Mat HF = Mat(width, height,CV_32FC1);

	for (double y = 0; y < height; y++)
	{
		for (double x = 0; x < width; x++)
		{
			double xn = GetXn(x, width);
			double yn = GetYn(y, height);

			double radius = sqrt(xn * xn + yn * yn);
			double ttmp2;

			if (radius <= x1) { ttmp2 = 0.0; }
			else if (radius >= x2) { ttmp2 = 1.0; }
			else // if ((radius > x1) && (radius < x2))
			{
				double ttmp = M_PI / 4.0 * (1.0 + (radius - x1) / (x2 - x1));
				ttmp2 = cos(M_PI / 2.0 * (log(ttmp * 2.0 / M_PI) / log(2.0)));
			}

			HF.at<float>(y, x) = ttmp2;
		}
	}
	
	Mat planes[] = { HF,HF };
	merge(planes, 2, HF);
	return HF;
}

Mat GenerateLF0_f(int width, int height, double x1, double x2)
{
	Mat LF = Mat(width, height, CV_32FC1);

	for (double y = 0; y < height; y++)
	{
		for (double x = 0; x < width; x++)
		{
			double xn = GetXn(x, width);
			double yn = GetYn(y, height);

			double radius = sqrt(xn * xn + yn * yn);

			if (radius < x1)
			{
				LF.at<float>(y, x) = 1.0;
			}

			if (radius > x2)
			{
				LF.at<float>(y, x) = 0.0;
			}

			if ((radius >= x1) && (radius <= x2))
			{
				double ttmp, ttmp2;
				ttmp = M_PI / 4.0 * (1.0 + (radius - x1) / (x2 - x1));
				ttmp2 = cos(M_PI / 2.0 * (log(ttmp * 4.0 / M_PI) / log(2.0)));

				LF.at<float>(y, x) = ttmp2;
			}
		}
	}
	Mat planes[] = { LF,LF };
	merge(planes, 2, LF);
	return LF;
}

Mat GenerateLF_f(int width, int height, double x1, double x2)
{
	Mat LF = Mat(width, height, CV_32FC1);

	for (double y = 0; y < height; y++)
	{
		for (double x = 0; x < width; x++)
		{
			double xn = GetXn(x, width);
			double yn = GetYn(y, height);

			double radius = sqrt(xn * xn + yn * yn);

			if (radius < x1)
			{
				LF.at<float>(y, x) = 1.0;
			}

			if (radius > x2)
			{
				LF.at<float>(y, x) = 0.0; 
			}

			if ((radius >= x1) && (radius <= x2))
			{
				double ttmp, ttmp2;
				ttmp = M_PI / 4.0 * (1.0 + (radius - x1) / (x2 - x1));
				ttmp2 = cos(M_PI / 2.0 * (log(ttmp * 4.0 / M_PI) / log(2.0)));

				LF.at<float>(y, x) = ttmp2; 
			}
		}
	}
	Mat planes[] = { LF,LF };
	merge(planes, 2, LF);
	return LF;
}


Mat GenerateDF_f(int width, int height, int k, int K)
{
	Mat DF_f = Mat(width, height, CV_32FC1);
	

	double alpha_K = pow(2.0,(double)K-1)*Factorial(K - 1.0) / sqrt(K * Factorial(2.0 * (K - 1.0)));

	
		for (double y = 0.0; y < height; y += 1.0)
		{
			for (double x = 0.0; x < width; x += 1.0)
			{				
				double xn = GetXn(x, width);
				double yn = GetYn(y, height);

				double theta = atan2(yn, xn) - (double)k * M_PI / K;
				double ttmp;
				if (abs(theta) < M_PI_2 || abs(theta) > 3 * M_PI_2)
				ttmp = alpha_K * pow((double)(cos(theta)), (double)(K - 1.0));
				else ttmp = -alpha_K * pow((double)(cos(theta)), (double)(K - 1.0));
				DF_f.at<float>(y, x) = ttmp;
			}
		}
		DF_f = DF_f/(float)(width * height);
		
		Mat planes[] = { DF_f,DF_f };
		merge(planes, 2, DF_f);
	return DF_f;
}