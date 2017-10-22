// BgEdgeDetect.h: interface for the BgEdgeDetect class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BGEDGEDETECT_H__EAEC8AA1_171F_4F85_933D_B9ED4E501D3D__INCLUDED_)
#define AFX_BGEDGEDETECT_H__EAEC8AA1_171F_4F85_933D_B9ED4E501D3D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "BgImage.h"
#include "BgEdgeList.h"

#define PI 3.1415926535
#define ZERO_TRESH 0.0000000001

// default values for edge detection
#define CONF_NMX 0.5
#define RANK_NMX 0.5
#define CONF_H 0.96
#define RANK_H 0.93
#define CONF_L 0.91
#define RANK_L 0.99
#define NMIN 5
#define KERNEL_SIZE2 2

#define HYST_LOW_CUT 0.0
#define MAX_CUSTT 30
#define MAX_FILTS 31
#define NO_ANGLES 361

#define ALF_TRESH PI/4

// 8邻域偏差
static const int gNb[8][2]=
{
   1, 0,
   1, 1,
   1,-1,
   0, 1,
   0,-1,
  -1, 0,
  -1, 1,
  -1,-1
};

// 8邻域角度
static const double gAlpha[8][2]=
{
     PI/2,   PI/2,
     PI/4,   PI/4,
   3*PI/4, 3*PI/4,
        0,   PI,
        0,   PI,
     PI/2,   PI/2,
   3*PI/4, 3*PI/4,
     PI/4,   PI/4
};

// main class, edge detection
class BgEdgeDetect
{
public:

   // main function for edge detection
   // cim input image
   // cel edge list (will be filled with pixels on edges)
   // nmxr, nmxc threshold for non-maxima-suppresion rank, confidence
   // rh, ch, threshold for hyst. high; rank, confidence
   // rl, cl, threshold for hyst. low; rank, confidence
   // nMin, min number of pixels on an edge
   // nmxType, hystTypeHigh, hystTypeLow, type of nmx curve, hyst. high curve, hyst low curve
   //  in (FC_ELLIPSE, FC_VERT_LINE, FC_HORIZ_LINE, FC_LINE, FC_SQUARE_BOX, FC_CUSTOM)

   void DoEdgeDetect(BgImage* cim, BgEdgeList* cel, double nmxr, double nmxc,
                     double rh, double ch, double rl, double cl,
                     int nMin, int nmxType, int hystTypeHigh, int hystTypeLow);

   // computes confidence map and rank information of sepcified image
   void ComputeEdgeInfo(BgImage*, float*, float*);
//   void ComputeConfidenceMap1(BgImage*, float*);
   // if have permanent data, call this function (compute only last two steps is same kernel size)
   void DoRecompute(BgEdgeList*, double, double, double, double, double, double, int, int, int, int);

   BgEdgeDetect(int filtDim);
   ~BgEdgeDetect();

   void SaveNmxValues();

   float EllipseEval(float, float);
   float EllipseComp(float, float, float, float);
   float LineEval(float, float);
   float LineComp(float, float, float, float);
   float VerticalLineEval(float, float);
   float VerticalLineComp(float, float, float, float);
   float HorizontalLineEval(float, float);
   float HorizontalLineComp(float, float, float, float);
   float SquareEval(float, float);
   float SquareComp(float, float, float, float);
   float CustomRegionEval(float, float);
   float CustomRegionComp(float, float, float, float);

   void SetCustomHigh(int*, int*, int, int, int);
   void SetCustomLow(int*, int*, int, int, int);
   void SetCustomHigh(double*, double*, int);
   void SetCustomLow(double*, double*, int);

   void IsGood(void);
   void GetPixels(int*, int*, int*, double, double, double, double);
   void GetNmxPixels(int*, int*, int*, double, double, double, double);
   
   double smofil_[MAX_FILTS];
   double diffil_[MAX_FILTS];
   double wdx_[MAX_FILTS*MAX_FILTS];
   double wdy_[MAX_FILTS*MAX_FILTS];
   double mN_[MAX_FILTS][MAX_FILTS];
   double mQ_[MAX_FILTS][MAX_FILTS];
   double* lookTable_[NO_ANGLES];

   int WW_;
   int WL_;
   float confTr_;
   float rankTr_;

   float* custx_;
   float* custy_;
   float* tcustx_;
   float* tcusty_;
   int ncust_;

   float* hcustx_;
   float* hcusty_;
   int nhcust_;
   float* lcustx_;
   float* lcusty_;
   int nlcust_;   

   int x_;
   int y_;
   float* permConf_;
   float* permRank_;
   float* permNmxRank_;
   float* permNmxConf_;
   bool havePerm_;

protected:

	// 创建角度蒙板
   void GenerateMaskAngle(double*, double); // 
   
   // 1.对图像进行高斯滤波，并提取梯度图像
   void CreateFilters(void);   // 创建滤波器
   void CreateLookTable(void); // 创建检查表
   void DeleteLookTable(void); // 删除检查表
   void GaussFilter(BgImage*, float*, double, int); // 高斯滤波
   void GaussDiffFilter(BgImage*, float*, float*, float*); // 高斯差分滤波
   
   // 2.进行新的非极大值抑制
   void Strength(float*, float*, float*); // 强度
   void NewNonMaxSupress(float*, float*, float*, float*, float*, float* ,
   float (BgEdgeDetect::*compf)(float, float, float, float)); // 非极大值抑制

   // 3.对置信度估计，并且进行新的磁滞效应和边界连接。对子空间进行估计
   void StrConfEstim(float*, float*, float*, float (BgEdgeDetect::*evalf)(float, float));
   void CompRanks(float*, float*); // 比较排名
   void NewHysteresisTr(float*, float*, BgEdgeList*, int, float*, float*); // 滞后效应
   void NewEdgeFollow(int, int); // 边界连接
   void SubspaceEstim(float*, float*, float*, float*);

   float* te_;
   float* tm_;
   double low_;
   float* tc_;
   float* tl_;
   int npt_;

   float* grx_;
   float* gry_;
   float* permGx_;
   float* permGy_; 
};

#endif // !defined(AFX_BGEDGEDETECT_H__EAEC8AA1_171F_4F85_933D_B9ED4E501D3D__INCLUDED_)
