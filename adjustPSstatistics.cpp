#include "stdafx.h"
#include "steerablepyramid.h"
#include "adjustPSstatistics.h"


Mat clip(Mat img, Point central, int width, int height)
{
	Mat cliped(img, Rect(central.x - (width-1) / 2, central.y - (height-1) / 2, width, height));
	return cliped;
}

Mat modAC(Mat noise, Mat Cy, float p)			//noise image is 1ch spacial image
{

	//calc AC of noise
	Mat Xf = MakeComplexMat(noise);
	dft(Xf, Xf);
	Mat Xf2;

	float size = Xf.cols * Xf.rows;
	int Nx = Xf.cols;
	int Ny = Xf.rows;
	int Nc = Cy.cols;							//size of AC

	vector<Mat> ch(2);
	split(Xf, ch);
	Mat power;
	power = ch[0].mul(ch[0]) + ch[1].mul(ch[1]);    //power spectrum
	Mat plane[] = { power,Mat::zeros(power.size(),CV_32FC1) };
	merge(plane, 2, Xf2);
	Mat Cx = f2s(Xf2);

	//unnormalize the previously normalized correlation
	Cy = Cy * size;

	//create the Truth AC at this p
	int cx = Nx / 2 ;
	int cy = Ny / 2 ;
	int Lc = (Nc - 1) / 2;
	Mat Cy0 = Cy;
	Mat cCx = DownSample(Cx, Nc, Nc);
	Cy = p * Cy + (1 - p) * cCx;
	Cy = Cy.clone();							// to make Cy continuous

	//build the matrix that performs the convolution Cy1=Tcx*Ch1
	Cx = DownSample(Cx, Nc+2*Lc, Nc+2*Lc);
	Mat Tcx;

	//build Tcx
	for (int y = Lc ; y < 3 * Lc+1; ++y)
	{
		for (int x = Lc ; x < 3 * Lc + 1; ++x)
		{
			Mat ccx = clip(Cx, Point(x, y), Nc, Nc).clone();
			Mat ccxi;
			flip(ccx, ccxi, -1);						//ccxi = ccxの上下左右反転(flip(A,B,-1))
			ccx = ccx + ccxi;
			ccx.at<float>(Lc + 1, Lc + 1) /= 2;			//真ん中だけ半分
			ccx = ccx.reshape(0, 1);							//1行に変形
			Tcx.push_back(ccx);						//Tcxの最後の行にccx追加
		}
	}

	//calc Ch1
	Mat Cy1 = Cy.reshape(0, 1);
	Cy1 = Cy1.t();
	Tcx = Tcx.clone();
	Mat inv = Tcx.inv(DECOMP_SVD);
	inv = inv.t();
	Mat Ch1 = inv * Cy1;
	Mat Ch = Ch1.t();
	Ch = Ch.reshape(1, Nc);						//正方形に戻す　reshape Ch1

	//conv X with Ch
	//auxの中心にChをコピー
	Mat aux = Mat::zeros(Xf.size(), CV_32FC1);
	Rect roi(Point(cx - Lc, cy - Lc), Size(Nc, Nc));
	Mat auxROI = aux(roi);

	Ch.copyTo(auxROI);

	ftshift(aux);
	aux = MakeComplexMat(aux);
	dft(aux, aux);
	split(aux, ch);									//ch[0]=spectrum
	Mat Chf = abs(ch[0]);
	sqrt(Chf, Chf);
	//1ch to 2ch
	Mat plane2[] = { Chf,Chf };
	merge(plane2, 2, Chf);			
	Mat Yf = Xf.mul(Chf);							//multiply in fourier domain
	Mat Y = f2s(Yf);

	return Y;
}