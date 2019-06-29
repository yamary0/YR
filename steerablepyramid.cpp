#include "stdafx.h"
#include "filtergenerator.h"
#include "steerablepyramid.h"

void visual(Mat img, String id)
{
	normalize(img, img, 0.0, 1.0, NORM_MINMAX);
	colormap(img, 7);
	imshow(id, img);
}
 Mat MakeComplexMat(Mat img)
{
	Mat ComplexImg;
	copyMakeBorder(img, img, 0, getOptimalDFTSize(img.rows) - img.rows, 0, getOptimalDFTSize(img.cols) - img.cols, BORDER_CONSTANT, Scalar::all(0.0));
	Mat planes[] = { img,Mat::zeros(img.size(),CV_32FC1) };
	merge(planes, 2, ComplexImg);
	return ComplexImg;
}

 void ftshift(Mat& ftimage)
 {
	 // regroup four quadrants of image so that the origin is in the center of the image
	 int cx = ftimage.cols / 2, cy = ftimage.rows / 2;
	 Mat q0(ftimage, Rect(0, 0, cx, cy)), q1(ftimage, Rect(cx, 0, cx, cy));
	 Mat q2(ftimage, Rect(0, cy, cx, cy)), q3(ftimage, Rect(cx, cy, cx, cy)), tmp;
	 q0.copyTo(tmp);
	 q3.copyTo(q0);
	 tmp.copyTo(q3);
	 q1.copyTo(tmp);
	 q2.copyTo(q1);
	 tmp.copyTo(q2);
 }

 Mat f2s(Mat ftimg)				//fourier (2ch) to spacial image(1ch)
 {
	 Mat temp = ftimg.clone();
	 ftshift(temp);
	 dft(temp, temp);

	 vector<Mat> real(2);
	 split(temp, real);

	 return real[0];
}

 Mat DownSample(Mat ftimg, int width, int hight)
 {
	 Mat DS(ftimg, Rect(ftimg.cols / 2 - width / 2, ftimg.rows / 2 - hight / 2, width, hight));
	 return DS;
 }

 void colormap(Mat &img,int n)
 {
	 n = (float)n;
	 img.forEach<float>([n](float &p, const int *position)->void {
		 
		 for (float i = 0; i < n; i++)
		 {
			 if (i/n<=p&&p <= (i + 1) / n)p = i/ n;
		 }
		 
		 });
 }

void SteerablePyramid::decompose() 
{
	

	//ft(original image) at first
	//real image(1ch) to complex image(2ch)

	src_f = MakeComplexMat(src_s);
	int width = src_f.cols;
	int height = src_f.rows;

	
	//dft
	dft(src_f, src_f);
	ftshift(src_f);

	//first Highpass Filtering
	Mat HF0_f = GenerateHF0_f(width, height, width / 4, width / 2);
	HR0_f = src_f.mul(HF0_f); 


	HF_f.emplace_back(HF0_f);


	//generate first Lowpass Filter

	Mat LF0_f = GenerateLF0_f(width, height, width / 4, width / 2);
	Mat LR0_f = src_f.mul(LF0_f); 
	

	LF_f.emplace_back(LF0_f);
	LR_f.emplace_back(LR0_f);
	
	//recursion part
	for (int n = 0; n < N; ++n)
	{
		//LR_f‚ÌÅŒã‚Ì—v‘f‚ðŽQÆ
		Mat LR = LR_f.back();					
		width = LR.cols;
		height = LR.rows;

		// generate HF
		Mat HF = GenerateHF_f(width, height, width / 8, width / 4);
		HF_f.emplace_back(HF);

		for (int k = 0; k < K; ++k)
		{
			// generate DF & BF
			Mat DF = GenerateDF_f(width, height, k, K);
			Mat BF = HF.mul(DF);								//ƒtƒBƒ‹ƒ^ƒŠƒ“ƒO‚ÉŽg‚¤BF
			Mat lf = LF_f.back();
			Mat rBF = BF.mul(lf);				//•Û‘¶—pBF(LF*HF*DF) LF‚ÌÅŒã‚Ì—v‘f‚ðŽQÆ

			OF_f.emplace_back(DF);
			BF_f.emplace_back(rBF);

			//apply BF
			Mat BR = LR.mul(BF);

			BR_f.emplace_back(BR);
		}
		// generate LF
		Mat LF = GenerateLF_f(width, height, width / 8, width / 4);

		// apply LF
		LR = LR.mul(LF);

		//downsampling
		
		Mat DS = DownSample(LR, LR.cols / 2, LR.rows / 2);
		Mat DSF = DownSample(LF, LF.cols / 2, LF.rows / 2);


		LR_f.emplace_back(DS);
		LF_f.emplace_back(DSF);

	}
	//fourier to spacial
	//1 piece
	HR0_s = f2s(HR0_f);
	
	//N+1 pieces
	for (int i = 0; i < N+1; ++i)
	{
		LR_s.emplace_back(f2s(LR_f[i]));
		LF_s.emplace_back(f2s(LF_f[i]));
		HF_s.emplace_back(f2s(HF_f[i]));
	}
	//N~K pieces
	for (int i = 0; i < N * K; i++)
	{
		BR_s.emplace_back(f2s(BR_f[i]));
		OF_s.emplace_back(f2s(OF_f[i]));
		BF_s.emplace_back(f2s(BF_f[i]));
	}
	
}

void SteerablePyramid::setMag()
{
	for (int i = 0; i < BR_s.size(); ++i)
	{
		Mat Mag = abs(BR_s[i] - mean(BR_s[i])[0]);
		Mag -= mean(Mag)[0];

		MB.emplace_back(Mag);
	}
}
void paste(Mat &dst, Mat src, int x, int y) 
{

	int width = src.cols;
	int height = src.rows;

	Mat roi_dst = dst(Rect(x, y, width, height));
	src.copyTo(roi_dst);
}
Mat SteerablePyramid::getALL(Domain mode)
{
	int col=0, row=0;
	for (int i = 0; i < N; ++i)
	{
		row += pow(0.5, i) * src_s.rows;
	}
	for (int j = 0; j < K; ++j)
	{

		col += src_s.cols;
	}

	Mat all(row, col, CV_32FC1, Scalar::all(0.5));
	if (mode == S)
	{
		row = 0;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < K; j++)
			{
				int num = i * K + j;
				normalize(BR_s[num], BR_s[num], 0.0, 1.0, NORM_MINMAX);
				paste(all, BR_s[num], j*BR_s[num].cols, row);
			}
			row += BR_s[i * K].rows;
		}
		return all;
	}
	else if (mode == F)
	{
		row = 0;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < K; j++)
			{
				int num = i * K + j;
				normalize(BR_f[num], BR_f[num], 0.0, 1.0, NORM_MINMAX);
				paste(all, BR_f[num], j * BR_f[num].cols, row);
			}
			row += BR_f[i * K].rows;
		}
		return all;
	}
	else if (mode == M)
	{
		row = 0;
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < K; j++)
			{
				int num = i * K + j;
				normalize(MB[num], MB[num], 0.0, 1.0, NORM_MINMAX);
				paste(all, MB[num], j * MB[num].cols, row);
			}
			row += MB[i * K].rows;
		}
		return all;
	}
	else
	{
		cout << "undefined domain" << endl;
	}
}
Mat SteerablePyramid::getBR(int n, int k, Domain mode)
{
	if (1 <= n && n <= N) 
	{
		if (1 <= k && k <= K) 
		{
			int num = (n - 1) * K + (k - 1);				//calculate index No.
			if (mode == S)
			{
				return BR_s[num];
			}
			else if (mode == F)
			{
				return BR_f[num];
			}
			else
			{
				cout << "undefined domain" << endl;
			}
		}
		else
		{
			cout << "k is between 1 ` " << K << endl;
		}
	}
	else
	{
		cout << "n is between 1 ` " << N << endl;
	}
	return Mat::zeros(100, 100, CV_32FC1);
}

Mat SteerablePyramid::getLR(int n, Domain mode)
{
	if (1 <= n && n <= N + 1)
	{
		int num = n - 1;
		if (mode == S)
		{
			return LR_s[n];
		}
		else if (mode == F)
		{
			return LR_f[n];
		}
		else
		{
			cout << "undefined domain" << endl;
		}
	}
	else
	{
		cout << "n is between 1 ` " << N+1 << endl;
	}
	return Mat::zeros(100, 100, CV_32FC1);
}

Mat SteerablePyramid::getMB(int n, int k)
{
	if (1 <= n && n <= N)
	{
		if (1 <= k && k <= K)
		{
			int num = (n-1) * K + (k-1) ;
			return MB[num];
		}
		else
		{
			cout << "k is between 1 ` " << K << endl;
		}
	}
	else
	{
		cout << "n is between 1 ` " << N << endl;
	}
	return Mat::zeros(100, 100, CV_32FC1);
}

Mat SteerablePyramid::getHF(int n,Domain mode)
{
	if (1 <= n && n <= N + 1)
	{
		int num = n;				//index No.
		if (mode == S)
		{
			return HF_s[num];
		}
		else if (mode == F)
		{
			vector<Mat> ch(2);
			split(HF_f[num], ch);
			return ch[0];
		}
		else
		{
			cout << "undefined domain" << endl;
		}
	}
	else
	{
		cout << "n is between 1 ` " << N + 1 << endl;
	}
	return Mat::zeros(100, 100, CV_32FC1);
}

Mat SteerablePyramid::getLF(int n,Domain mode)
{
	if (1 <= n && n <= N + 1)
	{
		int num = n;				//index No.
		if (mode == S)
		{
			return LF_s[num];
		}
		else if (mode == F)
		{
			vector<Mat> ch(2);
			split(LF_f[num], ch);
			return ch[0];
		}
		else
		{
			cout << "undefined domain" << endl;
		}
	}
	else
	{
		cout << "n is between 1 ` " << N + 1 << endl;
	}
	return Mat::zeros(100, 100, CV_32FC1);
}

Mat SteerablePyramid::getBF(int n, int k,Domain mode)
{
	if (1 <= n && n <= N)
	{
		if (1 <= k && k <= K)
		{
			int num = (n - 1) * K + (k - 1);				//index No.
			if (mode == S)
			{
				return BF_s[num];
			}
			else if (mode == F)
			{
				vector<Mat> ch(2);
				split(BF_f[num], ch);
				return ch[0];
			}
			else
			{
				cout << "undefined domain" << endl;
			}
		}
		else
		{
			cout << "k is between 1 ` " << K << endl;
		}
	}
	else
	{
		cout << "n is between 1 ` " << N << endl;
	}
	return Mat::zeros(100, 100, CV_32FC1);
}

Mat SteerablePyramid::getOF(int n, int k,Domain mode)
{
	if (1 <= n && n <= N)
	{
		if (1 <= k && k <= K)
		{
			int num = (n - 1) * K + (k - 1);				//index No.
			if (mode == S)
			{
				return OF_s[num];
			}
			else if (mode == F)
			{
				vector<Mat> ch(2);
				split(OF_f[num], ch);
				return ch[0];
			}
			else
			{
				cout << "undefined domain" << endl;
			}
		}
		else
		{
			cout << "k is between 1 ` " << K << endl;
		}
	}
	else
	{
		cout << "n is between 1 ` " << N << endl;
	}
	return Mat::zeros(100, 100, CV_32FC1);
}


//SteerablePyramidƒNƒ‰ƒX‚ÌƒRƒ“ƒXƒgƒ‰ƒNƒ^
SteerablePyramid::SteerablePyramid(Mat original, int Scale, int Orientation) 
{
	src_s = original;
	N = Scale;
	K = Orientation;
}