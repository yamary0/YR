#pragma once

Mat clip(Mat img, Point central, int width, int height);
Mat modAC(Mat noise, Mat AC_t, float p);