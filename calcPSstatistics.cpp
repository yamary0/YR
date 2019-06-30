#include "stdafx.h"
#include "steerablepyramid.h"
#include "calcPSstatistics.h"
//Expectation
double ex(Mat Statistic)
{
	double size = Statistic.cols * Statistic.rows;
	double expectation = sum(Statistic)[0] / size;
	return expectation;
}

//Marginal
void Marginal::setRange()
{
	minMaxLoc(img, &min, &max,NULL, NULL);
	if (min == 0 && max == 1)cout << "this image is nomalized" << endl;
}
void Marginal::setMeanVar()
{
	Scalar m, sd;
	meanStdDev(img, m, sd);
	mean = m[0];
	var = sd[0] * sd[0];
}
void Marginal::setSkewKurt()
{
	double size = img.cols * img.rows;
	Mat temp = (img - mean)/sqrt(var);
	Mat third, fourth;
	//skew
	pow(temp, 3, third);
	skew = cv::mean(third)[0];
	//kurt
	pow(temp, 4, fourth);
	kurt = cv::mean(fourth)[0];
}
double Marginal::getMin()
{
	return min;
}
double Marginal::getMax()
{
	return max;
}
double Marginal::getMean()
{
	return mean;
}
double Marginal::getVar()
{
	return var;
}
double Marginal::getSkew()
{
	return skew;
}
double Marginal::getKurt()
{
	return kurt;
}
Marginal::Marginal(Mat src, Mode id)
{
	img = src;
	if (id == ORG)
	{
		setRange();
		setMeanVar();
		setSkewKurt();
	}
	else if (id == HR)
	{
		setMeanVar();
	}
	else if (id == LR)
	{
		setMeanVar();
		setSkewKurt();
	}
	else
	{
		cout << "フッ、私は不必要なようだな..." << endl;
	}
}

double calcSpectral(Mat mg)
{
	Scalar m, sd;
	meanStdDev(mg, m,sd);
	return m[0];
}

//Linear Position & Energy Position
Mat calcAC(Mat img,int N) 
{
	img = MakeComplexMat(img);
	dft(img, img);
	ftshift(img);

	Mat ac;
	vector<Mat> ch(2);
	split(img, ch);
	Mat power;
	power = ch[0].mul(ch[0]) + ch[1].mul(ch[1]);    //power spectrum
	Mat plane[] = { power,Mat::zeros(power.size(),CV_32FC1) };
	merge(plane, 2, ac);
	ac = f2s(ac);
	ac /= ac.cols * ac.rows;
	
	//clip N*N central part
	ac = DownSample(ac, N, N);
	return ac;
}

//Linear Scale
vector<Mat> calcLS(Mat scl1, Mat scl2)		//scl1_s,scl2_s
{
	resize(scl2, scl2, Size(), 2.0, 2.0,INTER_NEAREST);

	scl1 = MakeComplexMat(scl1);
	scl2 = MakeComplexMat(scl2);

	dft(scl1, scl1);
	dft(scl2, scl2);
	ftshift(scl1);
	ftshift(scl2);
	vector<Mat> ch1, ch2;
	split(scl1, ch1);
	split(scl2, ch2);

	vector<Mat> im(3);
	im[0] =abs(ch1[0]);
	im[1] =abs(ch2[0]);
	im[2] =abs(ch2[1]);
	for (int i = 0; i < im.size(); i++)
	{
		visual(im[i], "fourier "+to_string(i));
	}

	Mat Rrp, Irp;
	Rrp = ch1[0].mul(ch2[0]);				//fine_r*curse_r
	Irp = ch1[0].mul(ch2[1]);				//fine_r*curse_i


	return {Rrp, Irp};
}

//Energy Orientation
Mat calcEO(Mat ori1, Mat ori2)		//input&output: spacial image
{
	Mat XC;
	XC = ori1.mul(ori2);
	return XC;
}

//Energy Scale
Mat calcES(Mat sca1, Mat sca2)		//sca1_s,sca2_s ,output:spacial image
{
	//double sca2_f phase
	resize(sca2, sca2, Size(), 2.0, 2.0,INTER_NEAREST);
	
	Mat XC;
	XC = sca1.mul(sca2);

	return XC;
}


/*
Mat RI2AP(Mat ftimg)					//RealImaginal to AmplitudePhase
{

	vector<Mat> ch(2);
	split(ftimg, ch);					//ch[0]:real ch[1]:imaginal
	Mat power, phase;

	//calc power spectrum
	power = ch[0].mul(ch[0]) + ch[1].mul(ch[1]); 
	sqrt(power, power);

	//calc phase
	phase = ch[1] / ch[0];
	phase.forEach<float>([](float& p, const int* position)-> void {
		p = atan2(p, 1);
		});

	//phase double
	phase *= 2;

}
*/