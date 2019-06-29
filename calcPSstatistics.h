#pragma once

//Šú‘Ò’l
double ex(Mat Statistic);

//Marginal
enum Mode {
	ORG,	//original
	HR,		//high pass residual
	LR		//low pass residual
};

class Marginal {
	protected:
		Mat img;
		double min;
		double max;
		double mean;
		double var;
		double skew;
		double kurt;
		void setRange();
		void setMeanVar();
		void setSkewKurt();
	public:
		double getMin();
		double getMax();
		double getMean();
		double getVar();
		double getSkew();
		double getKurt();
		Marginal(Mat src, Mode id);

};

double Spectral(Mat mg);					//Spectral ( mean of Magnitude )
Mat calcAC(Mat img,int N);					//Linear Position & Energy Position ( Auto-Correlation )
vector<Mat> calcLS(Mat scl1, Mat scl2);		//Linear Scale ( Rerative Phase )

Mat calcXC_Ori(Mat ori1, Mat ori2);			//Energy Orientation ( Cross-Correlation )
Mat calcXC_Sca(Mat sca1, Mat sca2);			//Energy Scale ( Cross-Correlation )

