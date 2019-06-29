#include "stdafx.h"
#include "steerablepyramid.h"
#include "calcPSstatistics.h"
#include "adjustPSstatistics.h"

constexpr int SCALE = 3;
constexpr int ORIENT = 4;
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
	cout << img.size() << endl;

	SteerablePyramid pyr(img, SCALE, ORIENT);
	Mat out;
	pyr.decompose();
	pyr.setMag();

	for (int i = 1; i < ORIENT + 1; i++) {
		Mat out = pyr.getBF(1, i, 'f');
		visual(out, to_string(i));
	}
	
	/*
	Mat all = pyr.getALL('m');
	normalize(all, out, 0.0, 1.0, NORM_MINMAX);
	out.convertTo(out, CV_8U, 255);
	imshow("all", all);
	//imwrite("all.png", out);

	Mat LR = pyr.getLR(1,'s');

	normalize(LR, out, 0.0, 1.0, NORM_MINMAX);
	out.convertTo(out, CV_8U,255);
	imwrite("LR.png", out);
	visual(LR, "LR");
	Mat AC = calcAC(LR, 7);
	normalize(AC, AC, 0.0, 1.0, NORM_MINMAX);
	Mat colored;
	AC.convertTo(colored, CV_8U,255);
	applyColorMap(colored, colored, COLORMAP_JET);
	resize(colored, colored, Size(), 30, 30, INTER_NEAREST);
	imshow("AC", colored);
	imwrite("AC.png", colored);



	Mat scale1 = pyr.getMag(1, 1);
	Mat scale2 = pyr.getMag(1, 2);
	Mat scale3 = pyr.getMag(1, 4);
	
	Mat XO = calcXC_Ori(scale1, scale2);
	cout <<"xc1 : "<< mean(XO)[0] << endl;
	Mat X1 = calcXC_Ori(scale1, scale3);
	cout << "xc2 : " << mean(X1)[0] << endl;
	normalize(XO, XO, 0.0, 1.0, NORM_MINMAX);
	XO.convertTo(colored, CV_8U, 255);
	applyColorMap(colored, colored, COLORMAP_JET);
	imshow("XO", colored);
	colored.convertTo(colored, CV_8U);
	imwrite("XO.png", colored);
	*/
	/*
	vector<Mat> fr = calcLS(scale1, scale2);
	
	for (int i = 0; i < fr.size(); i++) {
		
		visual(fr[i], "relative phase" + to_string(i));
	}
	*/


	waitKey();
	return 0;
}