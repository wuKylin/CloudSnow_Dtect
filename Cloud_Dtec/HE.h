/**
* @function EqualizeHist_Demo.cpp
* @brief Demo code for equalizeHist function
* @author OpenCV team
*/

#include "opencv2/highgui/highgui.hpp"  
#include "opencv2/imgproc/imgproc.hpp"  
#include <iostream>  
#include <stdio.h>  

using namespace cv;
using namespace std;

// Split RGB image into three channels and do HE...  
Mat equalizeChannelHist(const Mat & inputImage)
{
	if (inputImage.channels() >= 3)
	{
		vector<Mat> channels;
		split(inputImage, channels);

		Mat B, G, R;

		equalizeHist(channels[0], B);
		equalizeHist(channels[1], G);
		equalizeHist(channels[2], R);

		vector<Mat> combined;
		combined.push_back(B);
		combined.push_back(G);
		combined.push_back(R);

		Mat result;
		merge(combined, result);

		return result;
	}

	return Mat();
}

// Convert into Ycrcb, and do HE... 
Mat equalizeIntensityHist(const Mat & inputImage)
{
	if (inputImage.channels() >= 3)
	{
		Mat ycrcb;

		cvtColor(inputImage, ycrcb, COLOR_BGR2YCrCb);

		vector<Mat> channels;
		split(ycrcb, channels);

		equalizeHist(channels[0], channels[0]);

		Mat result;
		merge(channels, ycrcb);

		cvtColor(ycrcb, result, COLOR_YCrCb2BGR);

		return result;
	}

	return Mat();
}

void getGrayImageHistImage(const Mat & src, Mat & histImage)
{
	Mat hist;
	int histSize = 256;

	calcHist(&src, 1, 0, Mat(), hist, 1, &histSize, 0);
	normalize(hist, hist, 0, histImage.rows, NORM_MINMAX, CV_32F);

	histImage = Scalar::all(255);
	int binW = cvRound((double)histImage.cols / histSize);

	for (int i = 0; i < histSize; i++)
		rectangle(histImage, Point(i*binW, histImage.rows),
		Point((i + 1)*binW, histImage.rows - cvRound(hist.at<float>(i))),
		Scalar::all(0), -1, 8, 0);
}


//int main(int, char** argv)
//{
//	Mat src, dst;
//	Mat intensity_color_dst;
//	Mat channel_color_dst;
//
//	const char* source_gray_window = "Source Gray Image";
//	const char* equalized_gray_window = "Equalized Gray Image";
//	const char* source_color_window = "Source Color Image";
//	const char* equalized_intensity_color_window = "Equalized Intensity Color Image";
//	const char* equalized_channels_color_window = "Equalized Channels Color Image";
//
//	/// Load image  
//	src = imread(argv[1], 1);
//
//	if (src.empty())
//	{
//		cout << "Usage: ./Histogram_Demo <path_to_image>" << endl;
//		return -1;
//	}
//
//	/// color image intensity equalization  
//	{
//		intensity_color_dst = equalizeIntensityHist(src);
//
//		namedWindow(source_color_window, WINDOW_AUTOSIZE);
//		namedWindow(equalized_intensity_color_window, WINDOW_AUTOSIZE);
//
//		imshow(source_color_window, src);
//		imshow(equalized_intensity_color_window, intensity_color_dst);
//	}
//
//	/// color image each channel equalization  
//	{
//		channel_color_dst = equalizeChannelHist(src);
//		namedWindow(equalized_channels_color_window, WINDOW_AUTOSIZE);
//		imshow(equalized_channels_color_window, channel_color_dst);
//	}
//
//	/// gray image equalization  
//	{
//		cvtColor(src, src, COLOR_BGR2GRAY);
//		equalizeHist(src, dst);
//
//		namedWindow(source_gray_window, WINDOW_AUTOSIZE);
//		namedWindow(equalized_gray_window, WINDOW_AUTOSIZE);
//
//		imshow(source_gray_window, src);
//		imshow(equalized_gray_window, dst);
//
//		/// get source gray image Histogram  
//		Mat graySrc_histImage = Mat::ones(200, 260, CV_8U) * 255;
//		getGrayImageHistImage(src, graySrc_histImage);
//		imshow("source gray image histogram", graySrc_histImage);
//
//		/// get equalized gray image Histogram  
//		Mat grayDst_histImage = Mat::ones(200, 260, CV_8U) * 255;
//		getGrayImageHistImage(dst, grayDst_histImage);
//		imshow("Equalized gray image histogram", grayDst_histImage);
//	}
//
//	/// Wait until user exits the program  
//	waitKey(0);
//
//	return 0;
//}