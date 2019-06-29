#include "stdafx.h"
#include "steerablepyramid.h"
#include "calcPSstatistics.h"
#include "adjustPSstatistics.h"

constexpr int SCALE = 2;
constexpr int ORIENT = 2;

//when you get a component of pyramid, you should select
//scale is natural number upto SCALE
//orient is natural number upto ORIENT

int main() {
	
	Mat src(200, 200, CV_32FC1, Scalar::all(0));
	circle(src, Point(src.cols / 2, src.rows / 2), src.cols / 4, Scalar(1), FILLED);

	//imshow("original disk image", src);
	

	Mat img = imread("wood.png", IMREAD_GRAYSCALE);
	img.convertTo(img, CV_32FC1);
	resize(img, img, Size(),0.4, 0.4);
	imshow("original", img);

	SteerablePyramid pyr(img, SCALE, ORIENT);
	Mat out;
	pyr.decompose();
	pyr.setMag();

	//Marginal
	Marginal org(img, ORG);
	cout << "mean = " << org.getMean() << endl;
	cout << "var = " << org.getVar() << endl;
	cout << "sd = " << sqrt(org.getVar()) << endl;
	cout << "skew = " << org.getSkew() << endl;
	cout << "kurt = " << org.getKurt() << endl;
	cout << "min = " << org.getMin() << endl;
	cout << "max = " << org.getMax() << endl;

	waitKey();
	return 0;
}