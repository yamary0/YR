#pragma once


void visual(Mat img, String id);
Mat MakeComplexMat(Mat img);
void ftshift(Mat& ftimage);
Mat f2s(Mat ftimg);
Mat DownSample(Mat ftimg, int width, int hight);	//clip central part of img 
void paste(Mat& dst, Mat src, int x, int y);		//(x,y) is left-top of the Rect

void colormap(Mat &img,int n);

enum Domain
{
	S,			//spacial
	F,			//fourier
	M			//magnitude
};

class SteerablePyramid
{
	private:
		Mat src_s;							//分解する画像(CV_32FC1)
		Mat src_f;
		int N;							//スケール
		int K;							//方向

		Mat HR0_s, HR0_f;				//firstHighpassResidual _s：空間　_f：周波数(2ch)
		vector<Mat> LR_s, LR_f;			//LowpassResidual L_1,..,L_N+1 LRは最後のLRも含むためN+1
		vector<Mat> BR_s, BR_f;			//BandpassResidual B_11,.,B_1K,B_21,...,B_NK N=Scale,K=Orientation
		vector<Mat> MB;				//BR_s's amplitude(called magnitude)
		
		vector<Mat> HF_s, HF_f;			//HighpassFilter
		vector<Mat> LF_s, LF_f;			//LospassFilter
		vector<Mat> OF_s, OF_f;			//OrientedFilter
		vector<Mat> BF_s, BF_f;			//BandpassFilter BF = LF*HF*DF at the same scale

	public:
		void decompose();
		void setMag();
		Mat getALL(Domain mode);
		Mat getBR(int n, int k,Domain mode);		//get a subband at Scale:n Direction:k 
		Mat getLR(int n,Domain mode);				//get a subband at Scale:n
		Mat getMB(int n, int k);				//get a subband at Scale:n Direction:k
		Mat getHF(int n, Domain mode);
		Mat getLF(int n, Domain mode);
		Mat getBF(int n, int k,Domain mode);
		Mat getOF(int n, int k,Domain mode);
		SteerablePyramid(Mat original, int Scale, int Orientation);
};

