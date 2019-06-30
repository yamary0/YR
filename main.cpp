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

	imshow("original disk image", src);
	

	Mat img = imread("wood.png", IMREAD_GRAYSCALE);
	img.convertTo(img, CV_32FC1,1/255.0);
	resize(img, img, Size(),0.4, 0.4);
	cout << mean(src)[0] << endl << endl;;
	//imshow("original", img);

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


	waitKey();
	return 0;
}