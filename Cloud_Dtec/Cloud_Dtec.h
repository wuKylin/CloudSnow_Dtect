#include "gdal_priv.h"
#include "gdal.h"
#include"dirent.h"
#include <iostream>
#include <math.h>
#include<string>
#include<vector>
#include <iomanip>
#include<opencv2\opencv.hpp>
#include<opencv2\highgui\highgui.hpp>
#include<iostream>
#include"MeanShift.h"
#include"msImageProcessor.h"
#include"BgImage.h"
#include"BgEdge.h"
#include"BgEdgeDetect.h"
#include<fstream>
#include<Afxwin.h>
#include"HE.h"


using namespace::cv;
using namespace::std;

#define FeatureNum 45
#define FeatureNumCT 14
#define FeatureNumCL 17
#define FeatureNumCS 19
#define FeatureNumTL 23
#define FeatureNumTS 32
#define FeatureNumLS 14

#define thres_p_CS 0.2
Mat src, dst;
int height = 500; 
int width = 500;

CvSVMParams ctparams, clparams, csparams, tlparams, tsparams, lsparams;
CvSVM ctSVM,clSVM,csSVM,tlSVM,tsSVM,lsSVM;

ofstream PrecisionTxt;
/***************************** MeanShift���õ��ĵ��� **********************************/
class BgPointSet
{
public:
	int* x_;
	int* y_;
	int n_;
	//wxPen pen_;
	int type_; // 0 circle
	// 1 point
	BgPointSet()
	{
		x_ = NULL;
		y_ = NULL;
	};
	~BgPointSet();
	void CleanData()
	{
		if (x_ != NULL)	{ delete[]x_; }
		if (y_ != NULL)	{ delete[]y_; }
	};
	void SetPoints(int*, int*, int);
};

class segParam
{
public:
	msImageProcessor *iProc;
	//�����ͼ�����
	BgImage* cbgImage_;
	BgImage* filtImage_;
	BgImage* segmImage_;
	BgImage* whiteImage_;

	//point set used to draw boundaries
	BgPointSet* boundaries_; // ���ڴ洢�߽�
	BgPointSet* regionPts_;  // ���ڴ洢������������

	RegionList *regionList; // �洢�߽���Ϣ
	int *regionLables; // �洢����α�ʶ,��ÿ�����ض�Ӧ�ĸ������
	int *regionPC; // �洢������������������صĸ��� 
	float *regionData; // �洢ÿ������ε�����
	int regionCount; // ����εĸ���

	//������ر���
	BOOL *isCS;//�ָ������Ƿ�Ϊ��ѩ����

	//segmemtation parameters
	int		sigmaS, minRegion, kernelSize;
	float	sigmaR, aij, epsilon, *gradMap_, *confMap_, *weightMap_, *customMap_;
	bool		edgeParamsHaveChanged_, hasBoundaries_, buseWeightMap_;
	SpeedUpLevel	speedUpLevel_;
	float speedUpThreshold_;
	int hasImage_, hasFilter_, hasSegment_;

	int operation;//�жϲ�����ʽ
public:
	segParam()
	{
		iProc = NULL;
		cbgImage_ = NULL;
		filtImage_ = NULL;
		segmImage_ = NULL; 
		whiteImage_ = NULL;
		boundaries_ = NULL; 
		regionPts_ = NULL;
		regionList = NULL;
		regionLables = NULL;
		regionPC = NULL;
		regionData = NULL;
		isCS = NULL;
		gradMap_ = confMap_ = weightMap_ = customMap_ = NULL;
	};
};


/*****************************   ��������   **********************************/
/*****************************************************************************/
#define GLCM_DIS 2				 //�Ҷȹ��������ͳ�ƾ���
#define GLCM_CLASS 16			 //����Ҷȹ��������ͼ��Ҷ�ֵ�ȼ���
#define GLCM_ANGLE_HORIZATION 0  //ˮƽ
#define GLCM_ANGLE_VERTICAL   1  //��ֱ
#define GLCM_ANGLE_DIGONAL    2  //�Խ�
#define GRAY_HISTOGRAM_NUM    16 //�Ҷ�ֱ��ͼ������
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#define LTP_THRESH 5        //LTP�е�T��������������ĶԳƵ�����ҶȲ�����ж�
/*****************************************************************************/
typedef std::vector<CPoint> ptLineTable;


/**************  �ָ��ÿ�������Ӧ�������࣬��������SVM����   ******************/
struct SmpFeature
{
	int is_cloud;   //��¼���������ƻ��Ƿ��ƣ�1�ƣ�-1ѩ
	float NIR, R, G, B;      //NIR��ֵ
	int area_pnum;  //�������
	double fractalDimension;//��ά��
	double GLCM_features[4];
	double CGLCM_features[9];
	double Edge_features[4];
	double LTP_features[18];
	int** mgray_b, **mgray_g, **mgray_r, **mgray_nr;
	double **mgray_h, **mgray_i, **mgray_s,**SF;
	double v2nir, v2r, v2g, v2b;

	SmpFeature()
	{
		is_cloud = 0;
		NIR = R = G = B = 0.0;
		area_pnum = 0;
	};
};


void transRGBtoHIS(const int R,const int G,const int B,double &H,double &I,double &S)
{
	if (R == 0 && G == 0 && B == 0)
	{
		H = 0.0;
		I = 0.0;
		S = 0.0;
	}
	else{
		int tempR = R, tempG = G, tempB = B, temp;
		double tempH;
		if (tempR > tempG)
		{
			temp = tempR;
			tempR = tempG;
			tempG = temp;
		}
		if (tempG > tempB)
		{
			temp = tempG;
			tempG = tempB;
			tempB = temp;
		}
		if (tempR > tempG)
		{
			temp = tempR;
			tempR = tempG;
			tempG = temp;
		}
		I = 1.0 / 3.0*(R + G + B);
		S = 1.0 - 3.0 / (R + G + B)*tempR;
		tempH = acos(0.5*(tempR - tempG + tempR - tempB) / (((tempR - tempG)*(tempR - tempG) + (tempR - tempG)*(tempG - tempB))*((tempR - tempG)*(tempR - tempG) + (tempR - tempG)*(tempG - tempB))));
		if (S == 0)
		{
			H = 0;
		}
		else if (B > G)
		{
			H = 360 - tempH;
		}
		else{
			H = tempH;
		}
	}
}

void calSF(BYTE *pImgData1, int height, int width, double **SF)
{
	double **mgray_h, **mgray_i, **mgray_s;
	mgray_h = new double*[height];
	mgray_i = new double*[height];
	mgray_s = new double*[height];
	double maxI, minI, maxS, minS;
	int tempCount = 0;
	for (int i = 0; i < height; i++)
	{
		mgray_h[i] = new double[width];
		mgray_i[i] = new double[width];
		mgray_s[i] = new double[width];
		for (int j = 0; j < width; j++)
		{
			transRGBtoHIS(pImgData1[width*height*2+i*width + j], pImgData1[width*height + i*width + j], pImgData1[i*width + j],
				mgray_h[i][j], mgray_i[i][j], mgray_s[i][j]);
			if (pImgData1[width*height * 2 + i*width + j] != 0 || pImgData1[width*height + i*width + j] != 0 || pImgData1[i*width + j] != 0)
			{
				if (tempCount == 0)
				{
					maxI = minI = mgray_i[i][j];
					maxS = minS = mgray_s[i][j];
					tempCount++;
				}
				else{
					if (maxI <  mgray_i[i][j])maxI = mgray_i[i][j];
					if (minI >  mgray_i[i][j])minI = mgray_i[i][j];
					if (maxS <  mgray_s[i][j])maxS = mgray_s[i][j];
					if (minS >  mgray_s[i][j])minS = mgray_s[i][j];
				}
			}
		}
	}
	tempCount = 0;
	double minSF, maxSF;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			if (pImgData1[width*height * 2 + i*width + j] != 0 || pImgData1[width*height + i*width + j] != 0 || pImgData1[i*width + j] != 0)
			{
				mgray_i[i][j] = (mgray_i[i][j] - minI) / (maxI - minI);
				mgray_s[i][j] = (mgray_s[i][j] - minS) / (maxS - minS);
				SF[i][j] = (mgray_i[i][j] + 1.0) / (mgray_s[i][j] + 1.0);
				if (tempCount == 0)
				{
					maxSF = minSF = SF[i][j];
					tempCount++;
				}
				else{
					if (maxSF <  SF[i][j])maxSF = SF[i][j];
					if (minSF >  SF[i][j])minSF = SF[i][j];
				}
			}
		}
	}

	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			if (pImgData1[width*height * 2 + i*width + j] != 0 || pImgData1[width*height + i*width + j] != 0 || pImgData1[i*width + j] != 0)
			{
				SF[i][j] = (SF[i][j] - minSF) / (maxSF - minSF)*255.0;
			}
		}
	}

	for (int i = 0; i < height; i++)
	{
		delete[]mgray_h[i];
		delete[]mgray_i[i];
		delete[]mgray_s[i];
	}
	delete[]mgray_h;
	delete[]mgray_i;
	delete[]mgray_s;
}

void calFractalDimension(int **m_gray, int height, int width, double &rvalue){
	const int wf = 256;//�������ά���ı߳�
	const int c = 2;//����
	//��ͼ�����ݴ���mat
	Mat m_src(height, width, CV_8UC1);
	for (int i = 0; i < height; i++){
		for (int j = 0; j < width; j++){
			m_src.at<uchar>(i, j) = m_gray[i][j];
		}
	}
	//��mat��׼��Ϊwf*wf
	int h, w;
	if (height > width){
		h = wf;
		w = wf*width*1.0 / height;
	}
	else {
		w = wf;
		h = wf*height*1.0 / width;
	}

	Mat m_temp(h, w, CV_8UC1);
	Mat m_dst(wf, wf, CV_8UC1, cvScalar(0));
	Mat dstROI = m_dst(cv::Rect(0, 0, w, h));
	resize(m_src, m_temp, cv::Size(w, h));
	m_temp.copyTo(dstROI, m_temp);
	// 	imshow("m_temp", m_temp);
	// 	imshow("m_src", m_src);
	//	imshow("m_dst", m_dst);
	//  imwrite("m_dst.jpg", m_dst);
	//	imwrite("m_dst.jpg", m_dst);

	//**********************//
	//�������ά��
	//**********************//
	struct point_f{
		double x;
		double y;
	};
	vector<point_f> v_p;
	vector<point_f> v_p2;
	double aa;
	for (int i = 1; i < log(wf) / log(c); i++){
		int nn = pow(c, i);//һ�������Ӿ������
		int ww = wf / nn; //����߳�
		int count = 0;//ͳ�Ʒǿվ���ĸ���
		bool flag = false;
		point_f p;
		for (int j = 0; j< nn; j++){
			for (int k = 0; k < nn; k++){
				for (int m = 0; m < ww; m++){
					for (int n = 0; n < ww; n++){
						if (m_dst.at<uchar>(j*ww + m, k*ww + n) > 0){
							count++;
							flag = true;
							break;
						}
					}
					if (flag){
						flag = false;
						break;
					}
				}
			}
		}
		p.x = ww;
		p.y = count;
		v_p.push_back(p);
		p.x = -log(1.0 / nn);
		p.y = log(count);
		aa = log(count);
		v_p2.push_back(p);
	}

	double aver_x = 0.0;
	double aver_y = 0.0;
	double sum1 = 0.0;
	double sum2 = 0.0;
	for (int i = 0; i < v_p2.size(); i++){
		aver_x += v_p2[i].x;
		aver_y += v_p2[i].y;
		sum1 += v_p2[i].x*v_p2[i].y;
		sum2 += v_p2[i].x*v_p2[i].x;
	}
	aver_x /= v_p2.size();
	aver_y /= v_p2.size();

	rvalue = (sum1 - v_p2.size()*aver_x*aver_y) / (sum2 - v_p2.size()*aver_x*aver_x);

	rvalue -= 1;//��һ����0-1
}

double calContrastForTwoBand(float *m_gray, int heigh, int widt, int NbaseImg, int NcompImg)
{
	//ȡˮƽ����ֱ�������ԽǷ����ֵ��ƽ��ֵ
	int i, j;

	if (m_gray == NULL)
		return false;

	int * glcm1 = new int[GLCM_CLASS * GLCM_CLASS];
	int * glcm2 = new int[GLCM_CLASS * GLCM_CLASS];
	int * glcm3 = new int[GLCM_CLASS * GLCM_CLASS];
	int * glcm4 = new int[GLCM_CLASS * GLCM_CLASS];
	int * histImage = new int[widt * heigh * 3];

	if (NULL == glcm1 || NULL == glcm2 || NULL == glcm3 || NULL == glcm4 || NULL == histImage)
		return false;

	//�Ҷȵȼ���---��GLCM_CLASS���ȼ�
	for (i = 0; i < heigh*widt * 3; i++)
	{
		histImage[i] = (int)(m_gray[i] * GLCM_CLASS / 256);
	}

	//��ʼ����������
	for (i = 0; i < GLCM_CLASS; i++)
	{
		for (j = 0; j < GLCM_CLASS; j++)
		{
			glcm1[i * GLCM_CLASS + j] = 0;
			glcm2[i * GLCM_CLASS + j] = 0;
			glcm3[i * GLCM_CLASS + j] = 0;
			glcm4[i * GLCM_CLASS + j] = 0;
		}
	}
	

	//����Ҷȹ�������
	int k, l;
	//ˮƽ����
	for (i = 0; i <heigh; i++)
	{
		for (j = 0; j < widt; j++)
		{
			l = histImage[heigh*widt*NbaseImg + i * widt + j];
			if (j + GLCM_DIS >= 0 && j + GLCM_DIS < widt)
			{
				k = histImage[heigh*widt*NcompImg + i * widt + j + GLCM_DIS];
				glcm1[l * GLCM_CLASS + k]++;
			}
			if (j - GLCM_DIS >= 0 && j - GLCM_DIS < widt)
			{
				k = histImage[heigh*widt*NcompImg + i * widt + j - GLCM_DIS];
				glcm1[l * GLCM_CLASS + k]++;
			}
		}
	}
	//��ֱ����
	for (i = 0; i <heigh; i++)
	{
		for (j = 0; j < widt; j++)
		{
			l = histImage[heigh*widt*NbaseImg + i * widt + j];
			if (i + GLCM_DIS >= 0 && i + GLCM_DIS < heigh)
			{
				k = histImage[heigh*widt*NcompImg + (i + GLCM_DIS) * widt + j];
				glcm2[l * GLCM_CLASS + k]++;
			}
			if (i - GLCM_DIS >= 0 && i - GLCM_DIS < heigh)
			{
				k = histImage[heigh*widt*NcompImg + (i - GLCM_DIS) * widt + j];
				glcm2[l * GLCM_CLASS + k]++;
			}
		}
	}
	//�ԽǷ���
	for (i = 0; i <heigh; i++)
	{
		for (j = 0; j < widt; j++)
		{
			l = histImage[heigh*widt*NbaseImg + i * widt + j];

			if (j + GLCM_DIS >= 0 && j + GLCM_DIS < widt && i + GLCM_DIS >= 0 && i + GLCM_DIS < heigh)
			{
				k = histImage[heigh*widt*NcompImg + (i + GLCM_DIS) * widt + j + GLCM_DIS];
				glcm3[l * GLCM_CLASS + k]++;
			}
			if (j - GLCM_DIS >= 0 && j - GLCM_DIS < widt && i - GLCM_DIS >= 0 && i - GLCM_DIS < heigh)
			{
				k = histImage[heigh*widt*NcompImg + (i - GLCM_DIS) * widt + j - GLCM_DIS];
				glcm3[l * GLCM_CLASS + k]++;
			}
		}
	}
	//����һ���ԽǷ���
	for (i = 0; i <heigh; i++)
	{
		for (j = 0; j < widt; j++)
		{
			l = histImage[heigh*widt*NbaseImg + i * widt + j];

			if (j + GLCM_DIS >= 0 && j + GLCM_DIS < widt && i - GLCM_DIS >= 0 && i - GLCM_DIS < heigh)
			{
				k = histImage[heigh*widt*NcompImg + (i - GLCM_DIS) * widt + j + GLCM_DIS];
				glcm4[l * GLCM_CLASS + k]++;
			}
			if (j - GLCM_DIS >= 0 && j - GLCM_DIS < widt && i + GLCM_DIS >= 0 && i + GLCM_DIS < heigh)
			{
				k = histImage[heigh*widt*NcompImg + (i + GLCM_DIS) * widt + j - GLCM_DIS];
				glcm4[l * GLCM_CLASS + k]++;
			}
		}
	}

	//��������ֵ
	double sum = heigh * widt;
	double entropy = 0, energy = 0, contrast = 0, homogenity = 0;
	for (i = 0; i <GLCM_CLASS; i++)
	{
		for (j = 0; j < GLCM_CLASS; j++)
		{
			//�Աȶ�
			contrast += ((i - j) * (i - j) * glcm1[i * GLCM_CLASS + j]) / sum;
			contrast += ((i - j) * (i - j) * glcm2[i * GLCM_CLASS + j]) / sum;
			contrast += ((i - j) * (i - j) * glcm3[i * GLCM_CLASS + j]) / sum;
			contrast += ((i - j) * (i - j) * glcm4[i * GLCM_CLASS + j]) / sum;
			//����
			//energy += (glcm1[i * GLCM_CLASS + j] * glcm1[i * GLCM_CLASS + j]) / sum;
			//energy += (glcm2[i * GLCM_CLASS + j] * glcm2[i * GLCM_CLASS + j]) / sum;
			//energy += (glcm3[i * GLCM_CLASS + j] * glcm3[i * GLCM_CLASS + j]) / sum;
			//energy += (glcm4[i * GLCM_CLASS + j] * glcm4[i * GLCM_CLASS + j]) / sum;
			////��
			//if (glcm1[i * GLCM_CLASS + j] > 0)
			//	entropy -= (glcm1[i * GLCM_CLASS + j] * log10(double(glcm1[i * GLCM_CLASS + j]))) / sum;
			//if (glcm2[i * GLCM_CLASS + j] > 0)
			//	entropy -= (glcm2[i * GLCM_CLASS + j] * log10(double(glcm2[i * GLCM_CLASS + j]))) / sum;
			//if (glcm3[i * GLCM_CLASS + j] > 0)
			//	entropy -= (glcm3[i * GLCM_CLASS + j] * log10(double(glcm3[i * GLCM_CLASS + j]))) / sum;
			//if (glcm4[i * GLCM_CLASS + j] > 0)
			//	entropy -= (glcm4[i * GLCM_CLASS + j] * log10(double(glcm4[i * GLCM_CLASS + j]))) / sum;

			////һ����
			//homogenity += (1.0 / (1 + (i - j) * (i - j)) * glcm1[i * GLCM_CLASS + j]) / sum;
			//homogenity += (1.0 / (1 + (i - j) * (i - j)) * glcm2[i * GLCM_CLASS + j]) / sum;
			//homogenity += (1.0 / (1 + (i - j) * (i - j)) * glcm3[i * GLCM_CLASS + j]) / sum;
			//homogenity += (1.0 / (1 + (i - j) * (i - j)) * glcm4[i * GLCM_CLASS + j]) / sum;
		}
	}
	
	delete[] glcm1;
	delete[] glcm2;
	delete[] glcm3;
	delete[] glcm4;
	delete[] histImage;
	
	//��������ֵ
	return contrast / 4;
}

void calColorGLCM(float *m_gray, int height, int width, double *featureVector)
{
	//ȡˮƽ����ֱ�������ԽǷ����ֵ��ƽ��ֵ
	int i;

	if (m_gray == NULL)
		return;

	i = 0;

	featureVector[i++] = calContrastForTwoBand(m_gray, height, width, 0, 0);
	featureVector[i++] = calContrastForTwoBand(m_gray, height, width, 0, 1);
	featureVector[i++] = calContrastForTwoBand(m_gray, height, width, 0, 2);
	featureVector[i++] = calContrastForTwoBand(m_gray, height, width, 1, 0);
	featureVector[i++] = calContrastForTwoBand(m_gray, height, width, 1, 1);
	featureVector[i++] = calContrastForTwoBand(m_gray, height, width, 1, 2);
	featureVector[i++] = calContrastForTwoBand(m_gray, height, width, 2, 0);
	featureVector[i++] = calContrastForTwoBand(m_gray, height, width, 2, 1);
	featureVector[i++] = calContrastForTwoBand(m_gray, height, width, 2, 2);
}

//��ȡ�߽������
void GetBoundryPointIndex(int**m_gray,int height,int width,vector<int> &BPoint)
{
	for (int i = 0; i < height;i++)
	{
		for (int j = 0; j < width;j++)
		{
			if (m_gray[i][j] != 0)
			{
				if (i == 0 || i == height - 1 || j == 0 || j == width - 1)
				{
					BPoint.push_back(i*width + j);
					continue;
				}
				if (m_gray[i - 1][j] == 0 || m_gray[i - 1][j + 1] == 0 || m_gray[i][j + 1] == 0 || m_gray[i + 1][j + 1] == 0 ||
					m_gray[i + 1][j] == 0 || m_gray[i + 1][j - 1] == 0 || m_gray[i][j - 1] == 0 || m_gray[i - 1][j - 1] == 0)
				{
					//���ٿɼӣ�ȥ������Ҫ�߽�㡣 
					int count = 0;
					if (m_gray[i - 1][j] == 0)    count++;
					if (m_gray[i - 1][j + 1] == 0)count++;
					if (m_gray[i][j + 1] == 0)    count++;
					if (m_gray[i + 1][j + 1] == 0)count++;
					if (m_gray[i + 1][j] == 0)    count++;
					if (m_gray[i + 1][j - 1] == 0)count++;
					if (m_gray[i][j - 1] == 0)    count++;
					if (m_gray[i - 1][j - 1] == 0)count++;

					if (count == 1 && (m_gray[i - 1][j + 1] == 0 || m_gray[i + 1][j + 1] == 0 || m_gray[i + 1][j - 1] == 0 || m_gray[i - 1][j - 1] == 0))
						continue;

					BPoint.push_back(i*width + j);
					
				}						
			}
		}
	}
}

//��ȡ�߽�����������
void edgeTracing(Mat edgeMat, vector<int>& TPoint, size_t length)//length,num of edge points
{
	int height = edgeMat.rows;
	int width = edgeMat.cols;
	Mat edgeMatPlus1(height + 2, width + 2, CV_8U, Scalar(0));
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			edgeMatPlus1.at<uchar>(i + 1, j + 1) = edgeMat.at<uchar>(i, j);
		}
	}
	int xPos, yPos, xPosPre, yPosPre, xPosPlus, yPosPlus;

	//�ҵ���ʼ���ٵ�͵ڶ������ٵ�
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			xPosPlus = i + 1;
			yPosPlus = j + 1;
			if (edgeMatPlus1.at<uchar>(xPosPlus, yPosPlus) == 255)
			{
				TPoint.push_back(i*width + j);
				if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus - 1) == 255)
				{
					TPoint.push_back((i + 1)*width + j - 1);
					i = height;
					break;
				}
				else if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus) == 255)
				{
					TPoint.push_back((i + 1)*width + j);
					i = height;
					break;
				}
				else if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus + 1) == 255)
				{
					TPoint.push_back((i + 1)*width + j + 1);
					i = height;
					break;
				}
				else if (edgeMatPlus1.at<uchar>(xPosPlus, yPosPlus + 1) == 255)
				{
					TPoint.push_back(i*width + j + 1);
					i = height;
					break;
				}
				cerr << " Tracing failed. " << endl;
				break;
			}
		}

	}
	//����
	int TraceIndex = TPoint[1];
	int count = 1;
	do
	{
		xPos = TPoint[count] / width;
		yPos = TPoint[count] % width;
		xPosPlus = xPos + 1;
		yPosPlus = yPos + 1;
		xPosPre = TPoint[count - 1] / width;
		yPosPre = TPoint[count - 1] % width;
		if (xPos - xPosPre == 0 && yPos - yPosPre == 1)
		{

			if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus + 1) == 255)
			{
				TPoint.push_back((xPos + 1)*width + yPos + 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus, yPosPlus + 1) == 255)
			{
				TPoint.push_back(xPos*width + yPos + 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus - 1, yPosPlus + 1) == 255)
			{
				TPoint.push_back((xPos - 1)*width + yPos + 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus - 1, yPosPlus) == 255)
			{
				TPoint.push_back((xPos - 1)*width + yPos);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus - 1, yPosPlus - 1) == 255)
			{
				TPoint.push_back((xPos - 1)*width + yPos - 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
		}
		if (xPos - xPosPre == -1 && yPos - yPosPre == 1)
		{
			if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus + 1) == 255)
			{
				TPoint.push_back((xPos + 1)*width + yPos + 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus, yPosPlus + 1) == 255)
			{
				TPoint.push_back(xPos*width + yPos + 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus - 1, yPosPlus + 1) == 255)
			{
				TPoint.push_back((xPos - 1)*width + yPos + 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus - 1, yPosPlus) == 255)
			{
				TPoint.push_back((xPos - 1)*width + yPos);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus - 1, yPosPlus - 1) == 255)
			{
				TPoint.push_back((xPos - 1)*width + yPos - 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus, yPosPlus - 1) == 255)
			{
				TPoint.push_back(xPos*width + yPos - 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
		}
		if (xPos - xPosPre == -1 && yPos - yPosPre == 0)
		{

			if (edgeMatPlus1.at<uchar>(xPosPlus - 1, yPosPlus + 1) == 255)
			{
				TPoint.push_back((xPos - 1)*width + yPos + 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus - 1, yPosPlus) == 255)
			{
				TPoint.push_back((xPos - 1)*width + yPos);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus - 1, yPosPlus - 1) == 255)
			{
				TPoint.push_back((xPos - 1)*width + yPos - 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus, yPosPlus - 1) == 255)
			{
				TPoint.push_back(xPos*width + yPos - 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus - 1) == 255)
			{
				TPoint.push_back((xPos + 1)*width + yPos - 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
		}
		if (xPos - xPosPre == -1 && yPos - yPosPre == -1)
		{
			if (edgeMatPlus1.at<uchar>(xPosPlus - 1, yPosPlus + 1) == 255)
			{
				TPoint.push_back((xPos - 1)*width + yPos + 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus - 1, yPosPlus) == 255)
			{
				TPoint.push_back((xPos - 1)*width + yPos);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus - 1, yPosPlus - 1) == 255)
			{
				TPoint.push_back((xPos - 1)*width + yPos - 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus, yPosPlus - 1) == 255)
			{
				TPoint.push_back(xPos*width + yPos - 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus - 1) == 255)
			{
				TPoint.push_back((xPos + 1)*width + yPos - 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>((xPosPlus + 1), yPosPlus) == 255)
			{
				TPoint.push_back((xPos + 1)*width + yPos);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
		}
		if (xPos - xPosPre == 0 && yPos - yPosPre == -1)
		{

			if (edgeMatPlus1.at<uchar>(xPosPlus - 1, yPosPlus - 1) == 255)
			{
				TPoint.push_back((xPos - 1)*width + yPos - 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus, yPosPlus - 1) == 255)
			{
				TPoint.push_back(xPos*width + yPos - 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus - 1) == 255)
			{
				TPoint.push_back((xPos + 1)*width + yPos - 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus) == 255)
			{
				TPoint.push_back((xPos + 1)*width + yPos);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus + 1) == 255)
			{
				TPoint.push_back((xPos + 1)*width + yPos + 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
		}
		if (xPos - xPosPre == 1 && yPos - yPosPre == -1)
		{
			if (edgeMatPlus1.at<uchar>(xPosPlus - 1, yPosPlus - 1) == 255)
			{
				TPoint.push_back((xPos - 1)*width + yPos - 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus, yPosPlus - 1) == 255)
			{
				TPoint.push_back(xPos*width + yPos - 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus - 1) == 255)
			{
				TPoint.push_back((xPos + 1)*width + yPos - 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus) == 255)
			{
				TPoint.push_back((xPos + 1)*width + yPos);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus + 1) == 255)
			{
				TPoint.push_back((xPos + 1)*width + yPos + 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus, yPosPlus + 1) == 255)
			{
				TPoint.push_back(xPos*width + yPos + 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
		}
		if (xPos - xPosPre == 1 && yPos - yPosPre == 0)
		{

			if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus - 1) == 255)
			{
				TPoint.push_back((xPos + 1)*width + yPos - 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus) == 255)
			{
				TPoint.push_back((xPos + 1)*width + yPos);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus + 1) == 255)
			{
				TPoint.push_back((xPos + 1)*width + yPos + 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus, yPosPlus + 1) == 255)
			{
				TPoint.push_back(xPos*width + yPos + 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus - 1, yPosPlus + 1) == 255)
			{
				TPoint.push_back((xPos - 1)*width + yPos + 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
		}
		if (xPos - xPosPre == 1 && yPos - yPosPre == 1)
		{
			if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus - 1) == 255)
			{
				TPoint.push_back((xPos + 1)*width + yPos - 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus) == 255)
			{
				TPoint.push_back((xPos + 1)*width + yPos);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus + 1, yPosPlus + 1) == 255)
			{
				TPoint.push_back((xPos + 1)*width + yPos + 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus, yPosPlus + 1) == 255)
			{
				TPoint.push_back(xPos*width + yPos + 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus - 1, yPosPlus + 1) == 255)
			{
				TPoint.push_back((xPos - 1)*width + yPos + 1);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
			if (edgeMatPlus1.at<uchar>(xPosPlus - 1, yPosPlus) == 255)
			{
				TPoint.push_back((xPos - 1)*width + yPos);
				count++;
				TraceIndex = TPoint[count];
				continue;
			}
		}
		cerr << "Tracing 2 failed." << endl;

	} while (!(abs(TraceIndex / width - TPoint[0] / width) <= 1 && abs(TraceIndex % width - TPoint[0] % width) <= 1) || count <= length * 2 / 3);
	
}

//��������ֱ��ͼ
void calEdgeCurHis(int **m_gray, int height, int width,float *CurHistogram)
{
	Mat srcMatPlus(height + 2, width + 2, CV_8U, Scalar(0));
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			srcMatPlus.at<uchar>(i + 1, j + 1) = m_gray[i][j];
		}
	}
	BOOL isDone;
	int xPos, yPos;
	do{
		isDone = false;
		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				xPos = i + 1;
				yPos = j + 1;
				if (srcMatPlus.at<uchar>(xPos, yPos) > 0)
				{
					int count = 0;

					if (srcMatPlus.at<uchar>(xPos - 1, yPos) == 0)    count++;
					if (srcMatPlus.at<uchar>(xPos - 1, yPos + 1) == 0)count++;
					if (srcMatPlus.at<uchar>(xPos, yPos + 1) == 0)    count++;
					if (srcMatPlus.at<uchar>(xPos + 1, yPos + 1) == 0)count++;
					if (srcMatPlus.at<uchar>(xPos + 1, yPos) == 0)    count++;
					if (srcMatPlus.at<uchar>(xPos + 1, yPos - 1) == 0)count++;
					if (srcMatPlus.at<uchar>(xPos, yPos - 1) == 0)    count++;
					if (srcMatPlus.at<uchar>(xPos - 1, yPos - 1) == 0)count++;

					if (count == 7)
					{
						{
							srcMatPlus.at<uchar>(xPos, yPos) = 0;
							m_gray[i][j] = 0;
							isDone = true;
						}
					}
				}

			}
		}
	} while (isDone);

	vector<int> BPoint;
	GetBoundryPointIndex(m_gray, height, width, BPoint);

	int *Curvature = new int[BPoint.size()];
	Mat srcMatPlus2(height + 4, width + 4, CV_8U, Scalar(0));
	Mat srcMat(height, width, CV_8U, Scalar(0));
	Mat CurvatureLabel(height, width, CV_8U, Scalar(255));
	Mat edgeMat(height, width, CV_8U, Scalar(0));
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			if (m_gray[i][j]>0){
				srcMatPlus2.at<uchar>(i + 2, j + 2) = 1;
				srcMat.at<uchar>(i, j) = 255;
			}
		}
	}
	for (int i = 0; i < 10; i++){ CurHistogram[i] = 0.0; }
	for (size_t i = 0; i < BPoint.size(); i++)
	{
		xPos = (BPoint[i] / width) + 2;
		yPos = (BPoint[i] % width) + 2;
		Curvature[i] = 
			- srcMatPlus2.at<uchar>(xPos - 2, yPos - 2) - srcMatPlus2.at<uchar>(xPos - 2, yPos - 1) - srcMatPlus2.at<uchar>(xPos - 2, yPos) - srcMatPlus2.at<uchar>(xPos - 2, yPos + 1) - srcMatPlus2.at<uchar>(xPos - 2, yPos + 2)
			- srcMatPlus2.at<uchar>(xPos - 1, yPos - 2) - srcMatPlus2.at<uchar>(xPos - 1, yPos - 1) - srcMatPlus2.at<uchar>(xPos - 1, yPos) - srcMatPlus2.at<uchar>(xPos - 1, yPos + 1) - srcMatPlus2.at<uchar>(xPos - 1, yPos + 2)
			- srcMatPlus2.at<uchar>(xPos, yPos - 2) - srcMatPlus2.at<uchar>(xPos, yPos - 1) + 24 * srcMatPlus2.at<uchar>(xPos, yPos) - srcMatPlus2.at<uchar>(xPos, yPos + 1) - srcMatPlus2.at<uchar>(xPos, yPos + 2)
			- srcMatPlus2.at<uchar>(xPos + 1, yPos - 2) - srcMatPlus2.at<uchar>(xPos + 1, yPos - 1) - srcMatPlus2.at<uchar>(xPos + 1, yPos) - srcMatPlus2.at<uchar>(xPos + 1, yPos + 1) - srcMatPlus2.at<uchar>(xPos + 1, yPos + 2)
			- srcMatPlus2.at<uchar>(xPos + 2, yPos - 2) - srcMatPlus2.at<uchar>(xPos + 2, yPos - 1) - srcMatPlus2.at<uchar>(xPos + 2, yPos) - srcMatPlus2.at<uchar>(xPos + 2, yPos + 1) - srcMatPlus2.at<uchar>(xPos + 2, yPos + 2);
		Curvature[i] = abs(Curvature[i]-10);
		CurHistogram[Curvature[i]]++;
		CurvatureLabel.at<uchar>(BPoint[i] / width, BPoint[i] % width) = Curvature[i] / 9.0 * 255;
		edgeMat.at<uchar>(BPoint[i] / width, BPoint[i] % width) = 255;
	}
	for (int i = 0; i < 10; i++){ 
		CurHistogram[i] /= BPoint.size();
	}
	delete[] Curvature;
	//imwrite("Curvature.bmp", CurvatureLabel);
	//imwrite("srcmat.bmp", srcMat);
	//imwrite("edgemat.bmp", edgeMat);
}

int FFT(complex<double>*TD, complex<double>*FD, int r)//rΪlog2N������������  
{
	LONG    count;              // ����Ҷ�任����    
	int     i, j, k;              // ѭ������  
	int     bfsize, p;
	double  angle;              // �Ƕ�     
	complex<double> *W, *X1, *X2, *X;

	count = 1 << r;               // ���㸶��Ҷ�任����  

	// ������������洢��  
	W = new complex<double>[count / 2];
	X1 = new complex<double>[count];
	X2 = new complex<double>[count];

	// �����Ȩϵ��  
	for (i = 0; i < count / 2; i++)
	{
		angle = -i * PI * 2 / count;
		W[i] = complex<double>(cos(angle), sin(angle));
	}

	// ��ʱ���д��X1  
	memcpy(X1, TD, sizeof(complex<double>) * count);

	// ���õ����㷨���п��ٸ���Ҷ�任  
	for (k = 0; k < r; k++)//kΪ��������ļ���  
	{
		for (j = 0; j < 1 << k; j++)
		{
			bfsize = 1 << (r - k);//������������������  
			for (i = 0; i < bfsize / 2; i++)
			{
				p = j * bfsize;
				X2[i + p] = X1[i + p] + X1[i + p + bfsize / 2];
				X2[i + p + bfsize / 2] = (X1[i + p] - X1[i + p + bfsize / 2])
					* W[i * (1 << k)];
			}
		}
		X = X1;
		X1 = X2;
		X2 = X;
	}
	// ��������  
	for (j = 0; j < count; j++)
	{
		p = 0;
		for (i = 0; i < r; i++)
		{
			if (j&(1 << i))
			{
				p += 1 << (r - i - 1);
			}
		}
		FD[j] = X1[p];
	}

	delete W;
	delete X1;
	delete X2;
	return 0;
}

//����༶�ҳ�����Ҷ������
void calMCLFD(int **m_gray, const int height, const int width, double *aValueMCLFD)
{
	Mat srcMatPlus(height + 2, width + 2, CV_8U, Scalar(0));
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			srcMatPlus.at<uchar>(i + 1, j + 1) = m_gray[i][j];
		}
	}
	BOOL isDone;
	int xPos, yPos;
	do{
		isDone = false;
		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				xPos = i + 1;
				yPos = j + 1;
				if (srcMatPlus.at<uchar>(xPos, yPos) > 0)
				{
					int count = 0;

					if (srcMatPlus.at<uchar>(xPos - 1, yPos) == 0)    count++;
					if (srcMatPlus.at<uchar>(xPos - 1, yPos + 1) == 0)count++;
					if (srcMatPlus.at<uchar>(xPos, yPos + 1) == 0)    count++;
					if (srcMatPlus.at<uchar>(xPos + 1, yPos + 1) == 0)count++;
					if (srcMatPlus.at<uchar>(xPos + 1, yPos) == 0)    count++;
					if (srcMatPlus.at<uchar>(xPos + 1, yPos - 1) == 0)count++;
					if (srcMatPlus.at<uchar>(xPos, yPos - 1) == 0)    count++;
					if (srcMatPlus.at<uchar>(xPos - 1, yPos - 1) == 0)count++;

					if (count == 7 || count == 8)
					{
						{
							srcMatPlus.at<uchar>(xPos, yPos) = 0;
							m_gray[i][j] = 0;
							isDone = true;
						}
					}
				}

			}
		}
	} while (isDone);

	vector<int> BPoint;
	GetBoundryPointIndex(m_gray, height, width, BPoint);

	Mat edgeMat(height, width, CV_8U, Scalar(0));
	for (int i = 0; i < BPoint.size(); i++)
	{
		edgeMat.at<uchar>(BPoint[i] / width, BPoint[i] % width) = 255;
	}

	vector<int> TPoint;
	edgeTracing(edgeMat, TPoint, BPoint.size());

	//Mat traceMat(height, width, CV_8U, Scalar(0));
	//for (int i = 0; i < TPoint.size(); i++)
	//{
	//	traceMat.at<uchar>(TPoint[i] / width, TPoint[i] % width) = 255;
	//}
	//imwrite("edge.bmp", edgeMat);
	//imwrite("trace.bmp", traceMat);


	/*calculate MCLFD*/
	unsigned int uN;//�ȼ����ɢ�� N=2^t
	double *aLengthArc = new double[TPoint.size()];
	double dLength = 0.0;
	int xPrePos, yPrePos;
	for (int i = 0; i < TPoint.size(); i++)
	{
		if (i == 0)
		{
			xPos = TPoint[i] % width;
			yPos = TPoint[i] / width;
			xPrePos = TPoint[TPoint.size() - 1] % width;
			yPrePos = TPoint[TPoint.size() - 1] / width;
			aLengthArc[i] = sqrt((xPos - xPrePos)*(xPos - xPrePos) + (yPos - yPrePos)*(yPos - yPrePos));
		}
		else{
			xPos = TPoint[i] % width;
			yPos = TPoint[i] / width;
			xPrePos = TPoint[i - 1] % width;
			yPrePos = TPoint[i - 1] / width;
			aLengthArc[i] = sqrt((xPos - xPrePos)*(xPos - xPrePos) + (yPos - yPrePos)*(yPos - yPrePos));
		}
		dLength += aLengthArc[i];
	}
	int aIndex[16];
	aIndex[0] = 0;
	int j = 1;
	double dLengthTemp = 0.0;
	for (int i = 1; i < 16; i++)
	{
		for (; j < TPoint.size(); j++)
		{
			dLengthTemp += aLengthArc[j];
			if ((aLengthArc[j + 1] == 1.0&&fabs(dLengthTemp - i*dLength / 16.0)<=0.5) || (aLengthArc[j + 1] > 1.0&&fabs(dLengthTemp - i*dLength / 16.0) < sqrt(2)-0.5))
			{
				aIndex[i] = j;
				break;
			}
		}
	}
	int xNextPos, yNextPos;
	double aMCLength[16][8];
	for (int i = 0; i < 16; i++)
	{
		for (int j = 1; j <= 8; j++)
		{
			xPos = TPoint[aIndex[i]] % width;
			yPos = TPoint[aIndex[i]] / width;
			if (i + j < 16){
				xNextPos = TPoint[aIndex[i + j]] % width;
				yNextPos = TPoint[aIndex[i + j]] / width;
			}
			else{
				xNextPos = TPoint[aIndex[i + j - 16]] % width;
				yNextPos = TPoint[aIndex[i + j - 16]] / width;
			}
			aMCLength[i][j - 1] = sqrt((xPos - xNextPos)*(xPos - xNextPos) + (yPos - yNextPos)*(yPos - yNextPos)) / dLength;
		}
	}

	int count = 0;
	complex<double> Src[16];
	complex<double> dst[16];
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 16; j++)
		{
			Src[j] = complex<double>(aMCLength[j][i], 0);	
		}
		FFT(Src, dst, 4);
		for (int j = 0; j < 4; j++)
		{
			aValueMCLFD[count++] = sqrt(dst[j].real()*dst[j].real() + dst[j].imag()*dst[j].imag());
		}
	}
	delete[]aLengthArc;
}

int S_of_LTP(int dia_sub)
{
	if (dia_sub >= LTP_THRESH)
		return 2;
	if (dia_sub > -LTP_THRESH && dia_sub < LTP_THRESH)
		return 1;
	if (dia_sub <= -LTP_THRESH)
		return 0;
	return 0;
}

void Cal_IWCS_LTP(int **m_gray,const int height,const int width,double *iwcs_LTP)
{
	for (int i = 0; i < 18; i++)
	{
		iwcs_LTP[i] = 0.0;
	}
	int p_num = 0;
		for (int i = 0; i < height; i++)
		{
			for (int j = 0; j < width; j++)
			{
				if (i < 1 || i > height - 2 || j < 1 || j > width - 2|| m_gray[i][j]==0 )
					continue;
				p_num++;
				//���ֵÿ������9�ֿ��ܣ���������ֵ���Թ���һ��ʮ��ά����
				int iwcs_LTP_temp1 = S_of_LTP(m_gray[i - 1][j] - m_gray[i + 1][j]) * 1 +
					S_of_LTP(m_gray[i][j - 1] - m_gray[i][j + 1]) * 3;

				//�����費��Ҫ���Ը�ʲô����ֹ����أ�
				// 			double weight_1 =(double)((m_gray[i-1][j] - m_gray[i+1][j]) * 
				// 									  (m_gray[i-1][j] - m_gray[i+1][j]) +
				// 									  (m_gray[i][j-1] - m_gray[i][j+1]) *
				// 									  (m_gray[i][j-1] - m_gray[i][j+1])) / 130050.0;

				double weight_1 = (double)((m_gray[i - 1][j] - m_gray[i + 1][j]) *
					(m_gray[i - 1][j] - m_gray[i + 1][j]) +
					(m_gray[i][j - 1] - m_gray[i][j + 1]) *
					(m_gray[i][j - 1] - m_gray[i][j + 1]));

				int iwcs_LTP_temp2 = S_of_LTP(m_gray[i - 1][j - 1] - m_gray[i + 1][j + 1]) * 1 +
					S_of_LTP(m_gray[i + 1][j - 1] - m_gray[i - 1][j + 1]) * 3;
				double weight_2 = (double)((m_gray[i - 1][j - 1] - m_gray[i + 1][j + 1]) *
					(m_gray[i - 1][j - 1] - m_gray[i + 1][j + 1]) +
					(m_gray[i + 1][j - 1] - m_gray[i - 1][j + 1]) *
					(m_gray[i + 1][j - 1] - m_gray[i - 1][j + 1]));

				iwcs_LTP[iwcs_LTP_temp1] += weight_1;
				iwcs_LTP[iwcs_LTP_temp2 + 9] += weight_2;

				if (iwcs_LTP[iwcs_LTP_temp1] < 0 || iwcs_LTP[iwcs_LTP_temp2 + 9] < 0)
				{
					cerr << "LTPԪ������̫���ˣ�����ˣ�" << endl;
					return;
				}
			}
		}
		for (int j = 0; j<18; j++)
		{
			// ��һ��
			if (iwcs_LTP[j] > 0 && p_num != 0)   //svm_feature[i].iwcs_LTP[j]һ��Ϊ0��������
			{
				iwcs_LTP[j] = 1.0 / (1.0 + exp(log(iwcs_LTP[j] / p_num) / 10));
			}
		}
	
}