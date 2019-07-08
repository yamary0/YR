#include "stdafx.h"
#include "steerablepyramid.h"
#include "calcPSstatistics.h"
#include "adjustPSstatistics.h"

constexpr int SCALE = 2;
constexpr int ORIENT = 3;
constexpr int N = 7;

//when you get a component of pyramid, you should select
//scale is natural number upto SCALE
//orient is natural number upto ORIENT

int main() {
	

	Mat src(200, 200, CV_32FC1, Scalar::all(0));
	circle(src, Point(src.cols / 2, src.rows / 2), src.cols / 4, Scalar(1), FILLED);

	//imshow("original disk image", src);
	

	Mat img = imread("wood.png", IMREAD_GRAYSCALE);
	img.convertTo(img, CV_32FC1,1/255.0);
	resize(img, img, Size(),0.4, 0.4);
	cout << mean(src)[0] << endl << endl;;
	imshow("original", img);

	SteerablePyramid pyr(img, SCALE, ORIENT);
	Mat out;
	pyr.decompose();
	pyr.subSpacialMean();


	//calcMarginal
	Marginal org(img, ORG);
	Marginal hi(pyr.getHR(S), HR);
	vector<Marginal> lo;
	for (Mat x:pyr.LR_s)
	{
		Marginal low(x, LR);
		lo.emplace_back(low);
	}

	//set Magnitude Band and calc Spectral and sub mean
	pyr.setMag();
	vector<double> Spec;
	for (Mat x : pyr.MB)Spec.emplace_back(calcSpectral(x));
	pyr.subMagnitudeMean();

	//calc Linear Position
	vector<Mat> LP;
	for (Mat x : pyr.LR_s)LP.emplace_back(calcAC(x, N));

	//calc Energy Position
	vector<Mat> EP;
	for (Mat x : pyr.MB)EP.emplace_back(calcAC(x, N));

	//calc Linear Scale
	vector<Mat> LS_r, LS_i;
	for (int scale = 1; scale < SCALE; ++scale)
	{
		for (int fine = 1; fine < ORIENT+1; ++fine)
		{
			for (int curse = 1; curse < ORIENT+1; ++curse)
			{
				vector<Mat> temp = calcLS(pyr.getBR(scale, fine, S), pyr.getBR(scale + 1, curse, S));
				LS_r.emplace_back(temp[0]);
				LS_i.emplace_back(temp[1]);
			}
		}
	}

	//calc Energy Orientation


	//calc Energy Scale
	vector<Mat> ES;
	for (int scale = 1; scale < SCALE; ++scale)
	{
		for (int fine = 1; fine < ORIENT + 1; ++fine)
		{
			for (int curse = 1; curse < ORIENT + 1; ++curse)
			{
				Mat temp = calcES(pyr.getMB(scale, fine), pyr.getMB(scale + 1, curse));
				ES.emplace_back(temp);
			}
		}
	}


	//output
	cout << "MARGINAL SIZE = " << lo.size() << endl;
	cout << "SPECTRAL SIZE = " << Spec.size() << endl;
	cout << "LINEAR POSITION SIZE = " << LP.size() << endl;
	cout << "ENERGY POSITION SIZE = " << EP.size() << endl;
	cout << "LINEAR SCALE SIZE (REAL)= " << LS_r.size() << endl;
	cout << "LINEAR SCALE SIZE (IMAG)= " << LS_i.size() << endl;
	cout << "ENERGY SCALE SIZE = " << ES.size() << endl;

	waitKey();
	return 0;
}