// BgEdgeDetect.cpp: implementation of the BgEdgeDetect class.
//
//////////////////////////////////////////////////////////////////////


#include "BgEdgeDetect.h"
#include <math.h>
#include "BgImage.h"
#include "BgEdge.h"
#include "BgEdgeList.h"
#include "BgEdgeDetect.h"
#include "BgDefaults.h"
#include <stdio.h>



#define TOL_E 2.2e-8
//#define TOL_E 0.05
#define SIGN(x) (x<0)? -1:1;

static int my_sign(double val)
{
	if(val>TOL_E)
		return 1;
	if(val<-TOL_E)
		return -1;
	return 0;
}

extern double factorial(double);
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BgEdgeDetect::BgEdgeDetect(int filtDim)
{
   havePerm_ = false;
   WL_ = filtDim;
   WW_ = 2*WL_+1;
   nhcust_ = 0;
   nlcust_ = 0;
   tcustx_ = new float[MAX_CUSTT];
   tcusty_ = new float[MAX_CUSTT];
   CreateFilters(); // 创建过滤器
   CreateLookTable(); // lookTable是什么作用 
}

BgEdgeDetect::~BgEdgeDetect()
{
   if (havePerm_==true)
   {
      delete [] permGx_;
      delete [] permGy_;
      delete [] permConf_;
      delete [] permRank_;
      delete [] permNmxRank_;
      delete [] permNmxConf_;
   }
   if (nhcust_>0)
   {
      delete [] hcustx_;
      delete [] hcusty_;
   }
   if (nlcust_>0)
   {
      delete [] lcustx_;
      delete [] lcusty_;
   }
   delete [] tcustx_;
   delete [] tcusty_;
   DeleteLookTable();
}

void BgEdgeDetect::IsGood(void)
{
   /*if (havePerm_==true)
      //bgLog("good\n");
   else
      //bgLog("bad\n");*/
}


// 按照椭圆方程的写
float BgEdgeDetect::EllipseEval(float x, float y)
{
   return ((x*x)/(rankTr_*rankTr_)+(y*y)/(confTr_*confTr_)-1);
}

// 构建椭圆
float BgEdgeDetect::EllipseComp(float x0, float y0, float x, float y)
{
//   return (EllipseEval(x,y)-EllipseEval(x0,y0));
   return ((x*x-x0*x0)/(rankTr_*rankTr_)+(y*y-y0*y0)/(confTr_*confTr_));
}

// 
float BgEdgeDetect::LineEval(float x, float y)
{
   return (confTr_*x+rankTr_*y-confTr_*rankTr_);
}

float BgEdgeDetect::LineComp(float x0, float y0, float x, float y)
{
//   return (LineEval(x,y)-LineEval(x0,y0));
   return (confTr_*(x-x0)+rankTr_*(y-y0));
}

float BgEdgeDetect::VerticalLineEval(float x, float y)
{
   return (x-rankTr_);
}

float BgEdgeDetect::VerticalLineComp(float x0, float y0, float x, float y)
{
//   return (VerticalLineEval(x,y)-VerticalLineEval(x0,y0));
   return (x-x0);

}

float BgEdgeDetect::HorizontalLineEval(float x, float y)
{
   return(y-confTr_);
}

float BgEdgeDetect::HorizontalLineComp(float x0, float y0, float x, float y)
{
//   return (HorizontalLineEval(x,y)-HorizontalLineEval(x0,y0));
   return (y-y0);
}

float BgEdgeDetect::SquareEval(float x, float y)
{
   if ((x/rankTr_)>(y/confTr_))
      return(x-rankTr_);
   else
      return(y-confTr_);
}

float BgEdgeDetect::SquareComp(float x0, float y0, float x, float y)
{
//   return(SquareEval(x,y)-SquareEval(x0,y0));
   static float tret;
   tret = ((x/rankTr_)>(y/confTr_)) ? x-rankTr_ : y-confTr_;
   tret -= ((x0/rankTr_)>(y0/confTr_)) ? x0-rankTr_ : y0-confTr_;
   return tret;
}

// 
float BgEdgeDetect::CustomRegionEval(float r,float c)
{
   //evaluate user region function
   //returns -1 if inside +1 if outside
   
   if ((r+c)<=ZERO_TRESH)
      return -1;
   int i;
   int crossings=0;
   float x;
   
   //shift to origin
   for (i=0; i<ncust_; i++)
   {
      tcustx_[i]=custx_[i]-r;
      tcusty_[i]=custy_[i]-c;
   }
   
   for (i=0; i<(ncust_-1); i++)
   {
      if ( (tcusty_[i]  >0 && tcusty_[i+1]<=0) ||
           (tcusty_[i+1]>0 && tcusty_[i]  <=0) )
      {
         x = (tcustx_[i]*tcusty_[i+1]-tcustx_[i+1]*tcusty_[i])/(tcusty_[i+1]-tcusty_[i]);
         if (x>0)
            crossings++;
      }	   
   }		
   
   if ((crossings % 2) ==1)
      return -1;
   else
      return 1;
}

//  典型的区域比较
float BgEdgeDetect::CustomRegionComp(float r0, float c0, float r, float c)
{
   return 0;
}

// 产生角度蒙板，这个角度蒙板是如何产生的
void BgEdgeDetect::GenerateMaskAngle(double* a,double theta) {
   static int sflag;
   static int i,j,k;
   static double cval[4];
   static double corner[2][4];
   static double sinv,cosv;
   static double intrs[2][4];
   static int scor[4],nscor[4];
   static int sind,rowind,colind;
   static double cordi[2][4];
   static int lsigind,corin;
   static int sigind[4];
   static double diffin[2];
   static double comcoor;

   theta = theta*PI/180.0;// 转化为弧度制
   sinv = sin(theta); // 正弦
   cosv = cos(theta); // 余弦
   
   for (i=0; i<WW_*WW_; i++)
      a[i]=0;
   
   for (i=WL_; i>=-WL_; i--)
   {
      for(j=-WL_; j<=WL_; j++)
      {
         corner[0][0] = j-0.5;
         corner[0][1] = j+0.5;
         corner[0][2] = j+0.5;
         corner[0][3] = j-0.5;
         
         corner[1][0] = i+0.5;
         corner[1][1] = i+0.5;
         corner[1][2] = i-0.5;
         corner[1][3] = i-0.5;
         
         cval[0] = -sinv*corner[0][0]+cosv*corner[1][0];
         cval[1] = -sinv*corner[0][1]+cosv*corner[1][1];
         cval[2] = -sinv*corner[0][2]+cosv*corner[1][2];
         cval[3] = -sinv*corner[0][3]+cosv*corner[1][3];
         
         scor[0] = my_sign(cval[0]);
         scor[1] = my_sign(cval[1]);
         scor[2] = my_sign(cval[2]);
         scor[3] = my_sign(cval[3]);
         
         sind = 0;
         if (scor[0]!=0)
            nscor[sind++] = scor[0];
         if (scor[1]!=0)
            nscor[sind++] = scor[1];
         if (scor[2]!=0)
            nscor[sind++] = scor[2];
         if (scor[3]!=0)
            nscor[sind++] = scor[3];
         
         sflag = 0;
         for (k=1;k<sind;k++)
         {
            if (nscor[k]!=nscor[0])
               sflag++;
         }
         
         rowind = i+WL_;
         colind = j+WL_;
         
         if (sflag==0)
         {
            if (nscor[0]==1)
               a[colind+rowind*WW_] = 1.0;
            else
               a[colind+rowind*WW_] = 0.0;
         }
         
         if (sflag!=0)
         {
            for (k=0; k<4; k++)
               intrs[0][k] = intrs[1][k] = 0.0;
            
            if(scor[0]==0)
            {
               intrs[0][0] = corner[0][0];
               intrs[1][0] = corner[1][0];
            }
            if (scor[0]*scor[1]<0)
            {
               intrs[0][0] = corner[1][0]*cosv/sinv;
               intrs[1][0] = corner[1][0];
            }
            if (scor[1]==0)
            {
               intrs[0][1] = corner[0][1];
               intrs[1][1] = corner[1][1];
            }
            if (scor[1]*scor[2]<0)
            {
               intrs[0][1] = corner[0][1];
               intrs[1][1] = corner[0][1]*sinv/cosv;
            }
            if (scor[2]==0)
            {
               intrs[0][2] = corner[0][2];
               intrs[1][2] = corner[1][2];
            }
            if (scor[2]*scor[3]<0)
            {
               intrs[0][2] = corner[1][2]*cosv/sinv;
               intrs[1][2] = corner[1][2];
            }
            if (scor[3]==0)
            {
               intrs[0][3] = corner[0][3];
               intrs[1][3] = corner[1][3];
            }
            if (scor[3]*scor[0]<0)
            {
               intrs[0][3] = corner[0][3];
               intrs[1][3] = corner[0][3]*sinv/cosv;
            }
            
            corin = 0;
            if (fabs(intrs[0][0])>TOL_E || fabs(intrs[1][0])>TOL_E)
            {
               cordi[0][corin] = intrs[0][0];
               cordi[1][corin++] = intrs[1][0];
            }
            if (fabs(intrs[0][1])>TOL_E || fabs(intrs[1][1])>TOL_E)
            {
               cordi[0][corin] = intrs[0][1];
               cordi[1][corin++] = intrs[1][1];
            }
            if (fabs(intrs[0][2])>TOL_E || fabs(intrs[1][2])>TOL_E)
            {
               cordi[0][corin] = intrs[0][2];
               cordi[1][corin++] = intrs[1][2];
            }
            if (fabs(intrs[0][3])>TOL_E || fabs(intrs[1][3])>TOL_E)
            {
               cordi[0][corin] = intrs[0][3];
               cordi[1][corin++] = intrs[1][3];
            }
            
            lsigind=0;
			// 计算正数的个数
            if (scor[0]>0)
               sigind[lsigind++] = 0;
            if (scor[1]>0)
               sigind[lsigind++] = 1;
            if (scor[2]>0)
               sigind[lsigind++] = 2;
            if (scor[3]>0)
               sigind[lsigind++] = 3;
            
			
            if (lsigind==1)
            {
               a[colind+rowind*WW_] = 0.5*fabs(cordi[0][0]-cordi[0][1])*fabs(cordi[1][0]-cordi[1][1]);
            }
			
            if (lsigind==2)
            {
               diffin[0] = (int) fabs(cordi[0][0]-cordi[0][1]);
               diffin[1] = (int) fabs(cordi[1][0]-cordi[1][1]);
               if (diffin[0]==1)
               {
                  comcoor = corner[1][sigind[0]];
                  a[colind+rowind*WW_] = 0.5*(fabs(comcoor-cordi[1][0])+fabs(comcoor-cordi[1][1]));
               }
               if (diffin[1]==1)
               {
                  comcoor = corner[0][sigind[0]];
                  a[colind+rowind*WW_] = 0.5*(fabs(comcoor-cordi[0][0])+fabs(comcoor-cordi[0][1]));
               }
            }
		
            if(lsigind==3)
            {
				a[colind+rowind*WW_] = 1.0-0.5*fabs(cordi[0][0]-cordi[0][1])*fabs(cordi[1][0]-cordi[1][1]);
            }
         }
      }
   }
   
   // 对其进行规范化
   //A=A-mean(mean(A));
   comcoor = 0;
   for (i=0; i<WW_*WW_; i++)
      comcoor += a[i];
   comcoor /= WW_*WW_;
   for (i=0; i<WW_*WW_; i++)
      a[i] -= comcoor;
   
   //A=A/norm(A,'fro')
   comcoor = 0;
   for (i=0; i<WW_*WW_; i++)
      comcoor += a[i]*a[i];
   comcoor = sqrt(comcoor);
   for (i=0; i<WW_*WW_; i++)
      a[i] /= comcoor;
}

// 创建滤波器
void BgEdgeDetect::CreateFilters(void)
{
   int i,j;
   double w;
   for (i=-WL_; i<=WL_; i++)
   {
      w = pow((double)2,(-2*WL_))*factorial(2*WL_)/(factorial(WL_-i)*factorial(WL_+i));
      smofil_[i+WL_] = w;
      diffil_[i+WL_] = (2*i*w)/WL_;
   }
   for (j=0; j<WW_; j++)
   {
      for (i=0; i<WW_; i++)
      {
         wdy_[j+i*WW_] = wdx_[i+j*WW_] = smofil_[j]*diffil_[i];
      }
   }
   
   double norms = 0;
   double normd = 0;
   for (i=0; i<WW_; i++)
   {
      norms += smofil_[i]*smofil_[i];
      normd += diffil_[i]*diffil_[i];
   }
   
   for (j=0; j<WW_; j++)
   {
      for (i=0; i<WW_; i++)
      {
         mQ_[i][j] = (smofil_[j]*smofil_[i])/norms + (diffil_[j]*diffil_[i])/normd;
         mN_[i][j] = (i==j) ? 1-mQ_[i][j] : -mQ_[i][j];
      }
   }
}

void BgEdgeDetect::CreateLookTable()
{
   //bgLog("Creating angle lookup table\n");
   int i;
   for (i=-180; i<=180; i++)
   {
      lookTable_[i+180] = new double[WW_*WW_];
      GenerateMaskAngle(lookTable_[i+180], (double) i);
   }
}

void BgEdgeDetect::DeleteLookTable()
{
   int i;
   for (i=0; i<NO_ANGLES; i++)
   {
      delete [] lookTable_[i];
   }
}

// 获取像素
void BgEdgeDetect::GetPixels(int* nopix, int* pixx, int* pixy, double x1, double x2, double y1, double y2)
{
   double minx,maxx,miny,maxy;
   
   if (x1<x2)
   {
      minx = x1;
      maxx = x2;
   }
   else
   {
      minx = x2;
      maxx = x1;
   }
   
   if (y1<y2)
   {
      miny = y1;
      maxy = y2;
   }
   else
   {
      miny = y2;
      maxy = y1;
   }
   
   int i,j,npix;
   npix = 0;
   for (j=0; j<y_; j++)
   {
      for (i=0; i<x_; i++)
      {
         if (permRank_[i+j*x_]<maxx && permRank_[i+j*x_]>minx &&
             permConf_[i+j*x_]<maxy && permConf_[i+j*x_]>miny)
         {
            pixx[npix] = i;
            pixy[npix++] = j;
         }
      }
   }
   *nopix = npix;
}

// 获取非极大值像素
void BgEdgeDetect::GetNmxPixels(int* nopix, int* pixx, int* pixy, double x1, double x2, double y1, double y2)
{
   double minx,maxx,miny,maxy;
   if (x1<x2)
   {
      minx = x1;
      maxx = x2;
   }
   else
   {
      minx = x2;
      maxx = x1;
   }
   
   if (y1<y2)
   {
      miny = y1;
      maxy = y2;
   }
   else
   {
      miny = y2;
      maxy = y1;
   }
   
   int i,j,npix;
   npix = 0;
   for (j=0; j<y_; j++)
   {
      for (i=0; i<x_; i++)
      {
         if (permNmxRank_[i+j*x_]<maxx && permNmxRank_[i+j*x_]>minx &&
             permNmxConf_[i+j*x_]<maxy && permNmxConf_[i+j*x_]>miny)
         {
            pixx[npix] = i;
            pixy[npix++] = j;
	
	  }	 }
	  }
   *nopix = npix;
}

// 重新计算
void BgEdgeDetect::DoRecompute(BgEdgeList* cel, double nmxr, double nmxc,
                               double rh, double ch, double rl, double cl,
                               int nMin, int nmxType, int hystTypeHigh, int hystTypeLow)
{
   float *tr, *tc, *tdh, *tdl;
   //bgLog("Start edge detection...\n");
   tr = new float[x_*y_];
   tc = new float[x_*y_];
   tdh = new float[x_*y_];
   tdl = new float[x_*y_];
   
   //Nonmaximum supression
   //bgLog("...nonmaxima supression: ");

   float (BgEdgeDetect::*fcomp)(float,float,float,float);
   float (BgEdgeDetect::*feval)(float,float);
   
   switch(nmxType)
   {
   case FC_ELLIPSE:
      fcomp = &BgEdgeDetect::EllipseComp;
      feval = &BgEdgeDetect::EllipseEval;
      //bgLog("arc\n");
      break;
   case FC_VERT_LINE:
      fcomp = &BgEdgeDetect::VerticalLineComp;
      feval = &BgEdgeDetect::VerticalLineEval;
      //bgLog("vertical line\n");		
      break;
   case FC_HORIZ_LINE:
      fcomp = &BgEdgeDetect::HorizontalLineComp;
      feval = &BgEdgeDetect::HorizontalLineEval;
      //bgLog("horizontal line\n");		
      break;
   case FC_SQUARE_BOX:
      fcomp = &BgEdgeDetect::SquareComp;
      feval = &BgEdgeDetect::SquareEval;
      //bgLog("box\n");		
      break;
   case FC_LINE:
      fcomp = &BgEdgeDetect::LineComp;
      feval = &BgEdgeDetect::LineEval;		
      //bgLog("line\n");
      break;
  	default:
      //bgLog("Type not known\n");
      return;
   }
   confTr_ = (float) nmxc;
   rankTr_ = (float) nmxr;
   
   NewNonMaxSupress(permRank_,permConf_,permGx_,permGy_,permNmxRank_,permNmxConf_,fcomp);
   //bgLog("...hysteresis thresholding, high: ");
   
   switch(hystTypeHigh)
   {
   case FC_ELLIPSE:
     	fcomp = &BgEdgeDetect::EllipseComp;
      feval = &BgEdgeDetect::EllipseEval;
      //bgLog("arc");		
      break;
   case FC_VERT_LINE:
      fcomp = &BgEdgeDetect::VerticalLineComp;
      feval = &BgEdgeDetect::VerticalLineEval;
      //bgLog("vertical line");		
      break;
   case FC_HORIZ_LINE:
      fcomp = &BgEdgeDetect::HorizontalLineComp;
      feval = &BgEdgeDetect::HorizontalLineEval;
      //bgLog("horizontal line");		
      break;
   case FC_SQUARE_BOX:
      fcomp = &BgEdgeDetect::SquareComp;
      feval = &BgEdgeDetect::SquareEval;
      //bgLog("box");		
      break;
   case FC_LINE:
      fcomp = &BgEdgeDetect::LineComp;
      feval = &BgEdgeDetect::LineEval;		
      //bgLog("line");  		
      break;
  	case FC_CUSTOM:
      custx_ = hcustx_;
      custy_ = hcusty_;
      ncust_ = nhcust_;
      fcomp = &BgEdgeDetect::CustomRegionComp;
      feval = &BgEdgeDetect::CustomRegionEval;
      //bgLog("custom");
      break;
   }  
   confTr_ = (float) ch;
   rankTr_ = (float) rh;
   StrConfEstim(permNmxRank_, permNmxConf_, tdh, feval);

   //bgLog("  low: ");
   
   switch(hystTypeLow) {
   case FC_ELLIPSE:
     	fcomp = &BgEdgeDetect::EllipseComp;
      feval = &BgEdgeDetect::EllipseEval;
      //bgLog("arc\n");
      break;
   case FC_VERT_LINE:
      fcomp = &BgEdgeDetect::VerticalLineComp;
      feval = &BgEdgeDetect::VerticalLineEval;
      //bgLog("vertical line\n");		
      break;
   case FC_HORIZ_LINE:
      fcomp = &BgEdgeDetect::HorizontalLineComp;
      feval = &BgEdgeDetect::HorizontalLineEval;
      //bgLog("horizontal line\n");		
      break;
   case FC_SQUARE_BOX:
      fcomp = &BgEdgeDetect::SquareComp;
      feval = &BgEdgeDetect::SquareEval;
      //bgLog("box\n");		
      break;
   case FC_LINE:
      fcomp = &BgEdgeDetect::LineComp;
      feval = &BgEdgeDetect::LineEval;		
      //bgLog("line\n");  		
      break;
  	case FC_CUSTOM:
      custx_ = lcustx_;
      custy_ = lcusty_;
      ncust_ = nlcust_;
  	   fcomp = &BgEdgeDetect::CustomRegionComp;
  	   feval = &BgEdgeDetect::CustomRegionEval;
  	   //bgLog("custom\n");
  	   break;  		
   } 
   confTr_ = (float) cl;
   rankTr_ = (float) rl;
   StrConfEstim(permNmxRank_, permNmxConf_, tdl, feval);

  //hysteresis thresholding

   grx_ = permGx_;
   gry_ = permGy_;
   int minpt = nMin;
   NewHysteresisTr(tdh, tdl, cel, minpt, tr, tc);
   //bgLog("Done edge detection.\n");
   
   delete [] tdl;
   delete [] tdh;
   delete [] tr;
   delete [] tc;
}

// 保存非极大值抑制的值
void BgEdgeDetect::SaveNmxValues()
{
   FILE* fd;
   int i,j;
   
   fd = fopen("ranknmx.dat", "w");
   for (j=0; j<y_; j++)
   {
      for (i=0; i<x_; i++)
      {
         fprintf(fd, "%f ", *(permNmxRank_+j*x_+i));
      }
      fprintf(fd, "\n");
   }
   fclose(fd);
   
   fd=fopen("confnmx.dat", "w");
   for (j=0; j<y_; j++)
   {
      for (i=0; i<x_; i++)
      {
         fprintf(fd, "%f ", *(permNmxConf_+j*x_+i));
      }
      fprintf(fd, "\n");
   }
   fclose(fd);
}

//Computes confedence map and rank
//Pre : cim is an image
//Post: confidence map and rank has been computed for cim
//      and stored into confMap and rank respectively
// 计算边缘信息
void BgEdgeDetect::ComputeEdgeInfo(BgImage* cim, float* confMap, float *rank)
{
   x_ = cim->x_; // 图像的宽度
   y_ = cim->y_; // 图像的高度
   //bgLog("Computing confidence map...\n");
   
   // 计算图像的x，y方向的导数——通过高斯滤波差分
   float *pGx, *pGy, *pTemp;
   BgImage tcim(x_, y_);
   if (cim->colorIm_==1)
   {
	   tcim.SetImageFromRGB(cim->im_, x_, y_, false);
   } else
   {
	   tcim = *cim;
   }

   pGx = new float[x_*y_];
   pGy = new float[x_*y_];   
   pTemp = new float[x_*y_];
   
   // compute gradient images
   //bgLog("...smooth-differentiation filtering\n");
   GaussDiffFilter(&tcim, pGx, pGy, pTemp);   // 高斯差分滤波

   // compute confidences (subspace estimate)
   //bgLog("...subspace estimate\n");
   SubspaceEstim(pTemp, pGx, pGy, confMap);

   // compute edge strength from gradient image
   //bgLog("...edge strengths\n");
   Strength(pGx, pGy, pTemp); // 根据x，y方向的梯度求幅值
   
   // compute ranks of the strengths
   //bgLog("...computing ranks\n");
   CompRanks(pTemp, rank); // 计算排名
 
   //de-allocate memory
   delete [] pTemp;
   delete [] pGy;
   delete [] pGx;
}
/*
void BgEdgeDetect::ComputeConfidenceMap1(BgImage* cim, float* confMap)
{
	ComputeConfidenceMap(cim, confMap);
	BgEdgeList el;
	DoEdgeDetect(cim, &el, RANK_NMX, CONF_NMX, RANK_H, CONF_H, RANK_L, CONF_L,
		NMIN, FC_ELLIPSE, FC_SQUARE_BOX, FC_ELLIPSE);
   BgImage tempImage(cim->x_, cim->y_);
   el.SetBinImage(&tempImage);
   int i;
   for (i=0; i<(cim->x_*cim->y_); i++)
	   if (tempImage.im_[i] == 0)
		   confMap[i] = 0;
}
*/
void BgEdgeDetect::DoEdgeDetect(BgImage* cim, BgEdgeList* cel, double nmxr, double nmxc,
                                double rh, double ch, double rl, double cl,
                                int nMin, int nmxType, int hystTypeHigh, int hystTypeLow)
{
   x_ = cim->x_;
   y_ = cim->y_;
   //bgLog("Start edge detection...\n");   
   permGx_ = new float[x_*y_];
   permGy_ = new float[x_*y_];
   permConf_ = new float[x_*y_];
   permRank_ = new float[x_*y_];
   permNmxRank_ = new float[x_*y_];
   permNmxConf_ = new float[x_*y_];
   havePerm_ = true;
   float* tr;
   float* tc;
   float* tdh;
   float* tdl;
   
   tr = new float[x_*y_];
   tc = new float[x_*y_];
   tdh = new float[x_*y_];
   tdl = new float[x_*y_];
   
   // compute gradient images
   //bgLog("...smooth-differentiation filtering\n");
   GaussDiffFilter(cim, permGx_, permGy_, tr);   

   // compute confidences (subspace estimate)
   //bgLog("...subspace estimate\n");
   SubspaceEstim(tr, permGx_, permGy_, permConf_);
   
   // compute edge strength from gradient image
   //bgLog("...edge strengths\n");
   Strength(permGx_, permGy_, tr);
   
   // compute ranks of the strengths
   //bgLog("...computing ranks\n");
   CompRanks(tr, permRank_);
   
   // new nonmaxima supression
   //bgLog("...nonmaxima supression: ");

   // select appropriate function
   float (BgEdgeDetect::*fcomp)(float,float,float,float);
   float (BgEdgeDetect::*feval)(float,float);
   switch(nmxType)
   {
   case FC_ELLIPSE:
     	fcomp = &BgEdgeDetect::EllipseComp;
      feval = &BgEdgeDetect::EllipseEval;
      //bgLog("arc\n");
      break;
   case FC_VERT_LINE:
      fcomp = &BgEdgeDetect::VerticalLineComp;
      feval = &BgEdgeDetect::VerticalLineEval;
      //bgLog("vertical line\n");		
      break;
   case FC_HORIZ_LINE:
      fcomp = &BgEdgeDetect::HorizontalLineComp;
      feval = &BgEdgeDetect::HorizontalLineEval;
      //bgLog("horizontal line\n");		
      break;
   case FC_SQUARE_BOX:
      fcomp = &BgEdgeDetect::SquareComp;
      feval = &BgEdgeDetect::SquareEval;
      //bgLog("box\n");		
      break;
   case FC_LINE:
      fcomp = &BgEdgeDetect::LineComp;
      feval = &BgEdgeDetect::LineEval;		
      //bgLog("line\n");
      break;
  	default:
      //bgLog("Type not known\n");
      return;
   }

   confTr_ = (float) nmxc;
   rankTr_ = (float) nmxr;
   NewNonMaxSupress(permRank_, permConf_, permGx_, permGy_, permNmxRank_, permNmxConf_, fcomp);

   // new hysteresis thresholding
   //bgLog("...hysteresis thresholding, high: ");

   // select function, high curve
   switch(hystTypeHigh)
   {
   case FC_ELLIPSE:
     	fcomp = &BgEdgeDetect::EllipseComp;
      feval = &BgEdgeDetect::EllipseEval;
      //bgLog("arc");		
      break;
   case FC_VERT_LINE:
      fcomp = &BgEdgeDetect::VerticalLineComp;
      feval = &BgEdgeDetect::VerticalLineEval;
      //bgLog("vertical line");		
      break;
   case FC_HORIZ_LINE:
      fcomp = &BgEdgeDetect::HorizontalLineComp;
      feval = &BgEdgeDetect::HorizontalLineEval;
      //bgLog("horizontal line");		
      break;
   case FC_SQUARE_BOX:
      fcomp = &BgEdgeDetect::SquareComp;
      feval = &BgEdgeDetect::SquareEval;
      //bgLog("box");		
      break;
   case FC_LINE:
      fcomp = &BgEdgeDetect::LineComp;
      feval = &BgEdgeDetect::LineEval;		
      //bgLog("line");  		
      break;
  	case FC_CUSTOM:
      custx_ = hcustx_;
      custy_ = hcusty_;
      ncust_ = nhcust_;
      fcomp = &BgEdgeDetect::CustomRegionComp;
      feval = &BgEdgeDetect::CustomRegionEval;
      //bgLog("custom");
      break;
   }  

   confTr_ = (float) ch;
   rankTr_ = (float) rh;
   StrConfEstim(permNmxRank_, permNmxConf_, tdh, feval);

   //bgLog("  low: ");

   // select function, low curve
   switch(hystTypeLow)
   {
   case FC_ELLIPSE:
     	fcomp = &BgEdgeDetect::EllipseComp;
      feval = &BgEdgeDetect::EllipseEval;
      //bgLog("arc\n");
      break;
   case FC_VERT_LINE:
      fcomp = &BgEdgeDetect::VerticalLineComp;
      feval = &BgEdgeDetect::VerticalLineEval;
      //bgLog("vertical line\n");		
      break;
   case FC_HORIZ_LINE:
      fcomp = &BgEdgeDetect::HorizontalLineComp;
      feval = &BgEdgeDetect::HorizontalLineEval;
      //bgLog("horizontal line\n");		
      break;
   case FC_SQUARE_BOX:
      fcomp = &BgEdgeDetect::SquareComp;
      feval = &BgEdgeDetect::SquareEval;
      //bgLog("box\n");		
      break;
   case FC_LINE:
      fcomp = &BgEdgeDetect::LineComp;
      feval = &BgEdgeDetect::LineEval;		
      //bgLog("line\n");  		
      break;
  	case FC_CUSTOM:
      custx_ = lcustx_;
      custy_ = lcusty_;
      ncust_ = nlcust_;
  	   fcomp = &BgEdgeDetect::CustomRegionComp;
  	   feval = &BgEdgeDetect::CustomRegionEval;
  	   //bgLog("custom\n");
  	   break;  		
   } 
   confTr_ = (float) cl;
   rankTr_ = (float) rl;

   StrConfEstim(permNmxRank_, permNmxConf_, tdl, feval);

   grx_ = permGx_;
   gry_ = permGy_;
   
   NewHysteresisTr(tdh, tdl, cel, nMin, tr, tc);

   //bgLog("Done edge detection.\n");

   delete [] tdl;
   delete [] tdh;
   delete [] tr;
   delete [] tc;
}

// 按照x，y梯度计算置信度
// 置信度表示的是某处存在边缘的可信度，the presence of an edge. 是通过差分模板与标准化数据之间的相关系数求得，也就是他们夹角的余弦值
void BgEdgeDetect::SubspaceEstim(float* im, float* grx, float* gry, float* cee)
{
   // im original image
   // grx, gry gradient of image
   // cee confidence edge estimate
   
   float* itim;
   float* itgx;
   float* itgy;
   float* itcee;
   itim = im;
   itgx = grx;
   itgy = gry;
   itcee = cee;

   double* tae;
   double* ti;

   ti = new double[WW_*WW_];

   int i,j,l,c;
   double v1;
   double angleEdge;
   int WW2 = WW_*WW_;
   
   itim += WL_*x_;
   itgx += WL_*x_;
   itgy += WL_*x_;
   itcee += WL_*x_;

   // 按行进行，第WL_行开始
   for (j=WL_; j<(y_-WL_); j++)
   {
      for (i=0; i<WL_; i++)
         itcee[i] = 0;
      itim += WL_;
      itgx += WL_;
      itgy += WL_;
      itcee += WL_;
	  // 求置信度，
      for (i=WL_; i<(x_-WL_); i++, itim++, itgx++, itgy++, itcee++)
      {
         if ((fabs(*itgx)+fabs(*itgy))>TOL_E)
         {
            angleEdge = (-atan2(*itgx, *itgy))*180.0/PI; // 求边缘的方向角
            tae = lookTable_[(int) (angleEdge+180.49)]; 
            
            //A=A-mean(A)
            v1=0;
			// 取窗口，进行置信度计算
            for (l=0; l<WW_; l++)
            {
               for (c=0; c<WW_; c++)
               {
                  v1 += ti[l*WW_+c] = *(itim+(l-WL_)*x_+c-WL_);
               }
            }
            
			v1 /= WW2;
            for (l=0; l<WW2; l++)
               ti[l] -= v1;
            
            //A/norm(A,'fro')
            v1 = 0;
            for (l=0; l<WW2; l++)
               v1 += ti[l]*ti[l];
            v1 = sqrt(v1);
            for (l=0; l<WW2; l++)
               ti[l] /= v1; // 将向量进行归一化，标准化

            //global
            v1 = 0;
            for (l=0; l<WW2; l++)
               v1 += tae[l]*ti[l]; // 求角度模版与标准数据之间的相关性。
            v1 = fabs(v1);
            *itcee = (float) v1;
         }
         else
         {
            *itcee = 0;
         }
      }
      for (i=0; i<WL_; i++)
         itcee[i] = 0;
      itim += WL_;  // 图像信息
      itgx += WL_;  // x方向梯度求导
      itgy += WL_;  // y方向梯度求导
      itcee += WL_; // 置信度
   }
   WW2 = x_*y_;

   // 前后WL_的行置信度置为0
   for (j=0; j<(WL_*x_); j++)
   {
      cee[j] = 0;
      cee[WW2-j-1] = 0;
   }   

   
   delete [] ti;
}

// 高斯滤波器
void BgEdgeDetect::GaussFilter(BgImage* cim, float* fim, double sigma, int width)
{
   double* filter;
   unsigned char* im;
   double* tim;
   double sum = 0;
   double sum1 = 0;
   int i,ii,jj,j,k;

   im = cim->im_;
   if (width==-2)
   {
      for (i=0; i<x_*y_;i++)
         fim[i] = im[i];
      return;
   }

   if (width<3)
      width = (int) (1+2*ceil(2.5*sigma));
   int tail = width/2;
   width = 2*tail+1;

   //create kernel
   filter = new double[width];
   for (i=-tail; i<=tail; i++)
      sum += filter[i+tail] = exp(-i*i/(2*sigma*sigma));
   for (i=0; i<width; i++)
      filter[i] /= sum;

   //filter image
   im = cim->im_;
   tim = new double[x_*y_];
   for (j=0; j<y_; j++)
   {
	   // 进行一维卷积
      for (i=tail; i<(x_-tail); i++)
      {
         sum=0;
         for (k=-tail; k<=tail; k++)
            sum += im[j*x_+i+k]*filter[k+tail];
         tim[j*x_+i] = sum;
      }
      
      for (i=0; i<tail; i++)
      {
         tim[j*x_+i] = 0;
         tim[j*x_+x_-i-1] = 0;
         for (k=-tail; k<=tail; k++)
         {
            ii = (k+i)>=0 ? k+i : 0;
            tim[j*x_+i] += im[j*x_+ii]*filter[k+tail];
            ii = (x_-i-1+k)<x_ ? x_-i-1+k : x_-1;
            tim[j*x_+x_-i-1] += im[j*x_+ii]*filter[k+tail];
         }
      }
   }

   for (i=0; i<x_; i++)
   {
      for (j=tail; j<(y_-tail); j++)
      {
         sum=0;
         for (k=-tail; k<=tail; k++)
            sum += tim[(j+k)*x_+i]*filter[k+tail];
         fim[j*x_+i] = (float) (sum);
      }
      for (j=0; j<tail; j++)
      {
         sum = 0;
         sum1 = 0;
         for (k=-tail; k<=tail; k++)
         {
            jj = (k+j)>=0 ? k+j : 0;
            sum += tim[jj*x_+i]*filter[k+tail];
            jj = (y_-j-1+k)<y_ ? y_-j-1+k : y_-1;
            sum1 += tim[jj*x_+i]*filter[k+tail];
         }
         fim[j*x_+i] = (float) (sum);
         fim[(y_-j-1)*x_+i] = (float) (sum1);
      }
   }
   delete [] filter;
   delete [] tim;
}

void BgEdgeDetect::GaussDiffFilter(BgImage* cim, float* grx, float* gry, float* rezIm)
{
   
   double* sf; //smooth filter
   double* df; //diff filter
   unsigned char* im;
   
   double* tim;
   double sum = 0;
   double sum1 = 0;
   int i, j, k;
   
   //create kernels
   sf = smofil_; // 光滑滤波器
   df = diffil_; // 差分滤波器
   
   im = cim->im_;
   tim = new double[x_*y_];   //在这一句出错了，肯定是new的太多了
   for (i=0; i<x_*y_; i++)
   {
      grx[i] = gry[i] = 0;
      tim[i] = 0;
      rezIm[i] = im[i];
   }
   
   //filter image x
   //smooth on y
   for (i=0; i<x_; i++)
   {
      for (j=WL_; j<(y_-WL_); j++)
      {
         sum = 0;
         for (k=-WL_; k<=WL_; k++)
            sum += im[(j+k)*x_+i]*sf[k+WL_]; // 一维卷积运算（横向）
         tim[j*x_+i] = sum;
      }
   }
   
   //diff on x
   for (j=0; j<y_; j++)
   {
      for (i=WL_; i<(x_-WL_); i++)
      {
         sum = 0;
         for (k=-WL_; k<=WL_; k++)
            sum += tim[j*x_+i+k]*df[k+WL_]; // 纵向卷积运算
         grx[j*x_+i] = (float) (sum);
      }
   }

   //filter image y
   for (i=0; i<x_*y_;i++)
      tim[i] = 0;
   im = cim->im_;
   
   //smooth on x
   for (j=0; j<y_; j++)
   {
      for (i=WL_; i<(x_-WL_); i++)
      {
         sum = 0;
         for (k=-WL_; k<=WL_; k++)
            sum += im[j*x_+i+k]*sf[k+WL_];
         tim[j*x_+i] = sum;
      }
   }
  
   //diff on y
   for (i=0; i<x_; i++)
   {
      for (j=WL_; j<(y_-WL_); j++)
      {
         sum = 0;
         for (k=-WL_; k<=WL_; k++)
            sum += tim[(j+k)*x_+i]*df[k+WL_];
         gry[j*x_+i] = (float) (sum);
      }
   }
   delete [] tim;
}

// 比较排名，按照某点的边缘强度进行排名
void BgEdgeDetect::CompRanks(float* strength, float* ranks)
{
   int* index;
   float* ra;
   ra = new float[x_*y_];
   index = new int[x_*y_];
   int ii;
   
   for (ii=0; ii<x_*y_; ii++)
   {
      index[ii] = ii;
      ranks[ii] = 0;
      ra[ii] = strength[ii];
   }

   // 对排名进行堆排序
   //heap sort with ranks (from numerical recipies)
   unsigned long i, ir, j, l;
   unsigned long n;
   n = x_*y_;
   float rra;
   int irra;

   if (n<2)
      return;
   l = (n>>1)+1;
   ir = n;
   for (;;)
   {
      if (l>1)
      {
         rra = ra[--l-1];
         irra = index[l-1];
      }
      else
      {
         rra = ra[ir-1];
         irra = index[ir-1];
         ra[ir-1] = ra[1-1];
         index[ir-1] = index[1-1];
         if (--ir==1)
         {
            ra[1-1] = rra;
            index[1-1] = irra;
            break;
         }
      }
      i = l;
      j = l+l;


      while (j<=ir)
      {
         if (j<ir && ra[j-1]<ra[j+1-1])
            j++;
         if (rra<ra[j-1])
         {
            ra[i-1] = ra[j-1]; // 堆排序，将父亲结点赋予给孩子结点
            index[i-1] = index[j-1]; // 对标识进行相互交换
            i = j; // 将父亲结点赋予给孩子结点
            j <<= 1; // 对j进行乘以2，其作为父亲结点
         }
         else
            j = ir+1; // 
      }
      ra[i-1] = rra; // 对梯度强度进行排名
      index[i-1] = irra; // 记录
   }
   
   //setranks,设置排名
   irra = 1;
   for (ii=1; ii<x_*y_; ii++)
   {
      if (ra[ii]>ZERO_TRESH)
      {
         ranks[index[ii]] = (float) irra;
         if (ra[ii]>ra[ii-1])
            irra++;
      }
   }
   irra--;
   for (ii=0; ii<x_*y_; ii++) 
      ranks[ii] /= irra;

  delete [] index;
  delete [] ra;
}

// 进行置信度估计
void BgEdgeDetect::StrConfEstim(float* ranks, float* confidence, float* rezult,
                                float (BgEdgeDetect::*feval)(float,float))
{
   int i;
   for (i=0; i<x_*y_; i++)
   {
      rezult[i] = (this->*feval)(ranks[i], confidence[i]);
   }
}

// 某点处的边缘强度
void BgEdgeDetect::Strength(float* grx, float* gry, float* strength)
{
   int i,j;
   float* itgx;
   float* itgy;
   float* its;
   double val;
   itgx = grx;
   itgy = gry;
   its = strength;
   for (j=0; j<y_; j++)
      for(i=0; i<x_; i++)
      {
         val = sqrt(((double) (*itgx*(*(itgx++))))+((double) (*itgy*(*(itgy++)))));
         *(its++)=(float) (val);
      }
}

// 非极大值抑制
void BgEdgeDetect::NewNonMaxSupress(float* rank, float* conf, float* grx, float* gry, float* nmxRank, float* nmxConf,
                                    float (BgEdgeDetect::*feval)(float, float, float, float))
{
   int i,j;
   float* itr; // rank
   float* itc; // conf
   float* itgx;
   float* itgy;
   float* itnmxr; // 非极大值抑制rank
   float* itnmxc; // 非极大值抑制conf
   float alpha,r1,c1,r2,c2,lambda; // 梯度角度和幅值
   itr = rank;
   itc = conf;
   itgx = grx; // x方向的梯度
   itgy = gry; // y方向的梯度
   itnmxr = nmxRank; // 
   itnmxc = nmxConf; // 置信度confidence
   
   for (i=0; i<x_*y_; i++)
   {
      itnmxr[i] = itnmxc[i] = 0; 
   }

   // 四个边界
   for(i=0; i<x_; i++)
   {
      itr[i] = itc[i] = 0;
      itr[(y_-1)*x_+i] = itc[(y_-1)*x_+i] = 0;
   }
   for(j=0; j<y_; j++)
   {
      itr[j*x_] = itc[j*x_] = 0;
      itr[j*x_+x_-1] = itc[j*x_+x_-1] = 0;
   }

   // 3*3win
   for (j=0; j<y_; j++)
   {
      for (i=0; i<x_; i++, itr++, itc++, itgx++, itgy++, itnmxr++, itnmxc++)
      {
         if (*itr>0 && *itc>0)
         {
            alpha = (float) atan2(*itgy, *itgx); // 反三角求方向
            alpha = (alpha<0) ? alpha+(float)PI : alpha; // 转化为正值
			// 如果角度小于pi/4
            if (alpha<=PI/4)
            {
               lambda = (float) tan(alpha);
               r1 = (1-lambda)*(*(itr+1))+lambda*(*(itr+x_+1));
               c1 = (1-lambda)*(*(itc+1))+lambda*(*(itc+x_+1));
               r2 = (1-lambda)*(*(itr-1))+lambda*(*(itr-x_-1));
               c2 = (1-lambda)*(*(itc-1))+lambda*(*(itc-x_-1));
               if ((this->*feval)(*itr, *itc, r1, c1)<0 && (this->*feval)(*itr, *itc, r2, c2)<=0)
               {
                  *itnmxr = *itr;
                  *itnmxc = *itc;
               }
            }
            else if (alpha<=PI/2)
            {
               lambda = (float) tan(PI/2-alpha);
               r1 = (1-lambda)*(*(itr+x_))+lambda*(*(itr+x_+1));
               c1 = (1-lambda)*(*(itc+x_))+lambda*(*(itc+x_+1));
               r2 = (1-lambda)*(*(itr-x_))+lambda*(*(itr-x_-1));
               c2 = (1-lambda)*(*(itc-x_))+lambda*(*(itc-x_-1));
               if ((this->*feval)(*itr, *itc, r1, c1)<0 && (this->*feval)(*itr, *itc, r2, c2)<=0)
               {
                  *itnmxr = *itr;
                  *itnmxc = *itc;
               }
               
            }
            else if(alpha<=3*PI/4)
            {
               lambda = (float) tan(alpha-PI/2);
               r1 = (1-lambda)*(*(itr+x_))+lambda*(*(itr+x_-1));
               c1 = (1-lambda)*(*(itc+x_))+lambda*(*(itc+x_-1));
               r2 = (1-lambda)*(*(itr-x_))+lambda*(*(itr-x_+1));
               c2 = (1-lambda)*(*(itc-x_))+lambda*(*(itc-x_+1));
               if ((this->*feval)(*itr, *itc, r1, c1)<0 && (this->*feval)(*itr, *itc, r2, c2)<=0)
               {
                  *itnmxr = *itr;
                  *itnmxc = *itc;
               }
            }
            else
            {
               lambda = (float) tan(PI-alpha);
               r1 = (1-lambda)*(*(itr-1))+lambda*(*(itr+x_-1));
               c1 = (1-lambda)*(*(itc-1))+lambda*(*(itc+x_-1));
               r2 = (1-lambda)*(*(itr+1))+lambda*(*(itr-x_+1));
               c2 = (1-lambda)*(*(itc+1))+lambda*(*(itc-x_+1));
               if ((this->*feval)(*itr, *itc, r1, c1)<0 && (this->*feval)(*itr, *itc, r2, c2)<=0)
               {
                  *itnmxr = *itr;
                  *itnmxc = *itc;
               }
            }
         }
      }
   }
}

// 磁滞效应
void BgEdgeDetect::NewHysteresisTr(float* edge, float* low, BgEdgeList* cel, int nMin, float* mark, float* coord)
{
   float* tm;
   float* te;
   float* tl;
   int i,j,n;
   
   n=0;
   for (i=0, tm=mark; i<x_*y_; i++,tm++)
      *tm=0;
   
   te_ = te = edge;
   tm_ = tm = mark;
   tl_ = tl = low;
   
   for (j=0; j<y_; j++)
   {
      for (i=0; i<x_; i++, tm++, te++)
      {
		  // 将非极大值抑制的值并且大于高阈值的值的点开始 进行边界连接
         if ((*tm==0) && ((*te)>HYST_LOW_CUT))
         {
            //found an edge start
            npt_ = 0;
            *tm = 1;
            tc_ = coord;
            NewEdgeFollow(i, j);
            //store the edge
            if (npt_>=nMin) cel->AddEdge(coord, npt_);
         }
      }
   }
}

// 进行边缘连接
void BgEdgeDetect::NewEdgeFollow(int ii,int jj)
{
   int i;
   int iin, jjn;
   for (i=0; i<8; i++)
   {
      iin = ii+gNb[i][0];
      jjn = jj+gNb[i][1];
      if ((tm_[jjn*x_+iin]==0) && ((tl_[jjn*x_+iin])>0))
      {
         tm_[jjn*x_+iin] = 1;
         NewEdgeFollow(iin, jjn);
      }
   }
   *(tc_++) = (float) ii;
   *(tc_++) = (float) jj;
   npt_++;
}

// 设置高阈值
void BgEdgeDetect::SetCustomHigh(int* x, int* y, int n, int sx, int sy)
{
   if (nhcust_>0)
   {
      delete [] hcustx_;
      delete [] hcusty_;
   }
   nhcust_ = 0;
   hcustx_ = hcusty_ = 0;   
   nhcust_ = n+2;
   hcustx_ = new float[nhcust_];
   hcusty_ = new float[nhcust_];
   
   int idx,i;
   idx = 0;
   hcustx_[idx] = 0;
   hcusty_[idx++] = 0;
   for (i=0; i<n; i++)
   {
      hcustx_[idx] = ((float) x[i])/sx;
      hcusty_[idx++] = (float)(1.0-((float) y[i])/sy);
   }   
   hcustx_[idx] = 0;
   hcusty_[idx++] = 0;
   
   /*bgLog(" hyst high custom x: ");
   for (i=0; i<=n; i++)
      bgLog(" %f", hcustx_[i]);
   bgLog("\n");
   bgLog(" hist high custom y: ");
   for (i=0; i<=n; i++)
      bgLog(" %f", hcusty_[i]);
   bgLog("\n");*/   
}

// 设置高阈值
void BgEdgeDetect::SetCustomHigh(double* x, double* y, int n)
{
   if (nhcust_>0)
   {
      delete [] hcustx_;
      delete [] hcusty_;
   }
   nhcust_ = 0;
   hcustx_ = hcusty_ = 0;
   nhcust_ = n+2;
   hcustx_ = new float[nhcust_];
   hcusty_ = new float[nhcust_];
   
   int idx,i;
   idx = 0;
   hcustx_[idx] = 0;
   hcusty_[idx++] = 0;
   for (i=0; i<n; i++)
   {
      hcustx_[idx] = (float) x[i];
      hcusty_[idx++] = (float) y[i];
   }   
   hcustx_[idx] = 0;
   hcusty_[idx++] = 0;
   
   /*bgLog(" hyst high custom x: ");
   for (i=0; i<=n; i++)
      bgLog(" %f", hcustx_[i]);
   bgLog("\n");
   bgLog(" hist high custom y: ");
   for (i=0; i<=n; i++)
      bgLog(" %f", hcusty_[i]);
   bgLog("\n");*/   
}

// 设置低阈值
void BgEdgeDetect::SetCustomLow(int* x, int* y, int n, int sx, int sy)
{
   if(nlcust_>0)
   {
      delete [] lcustx_;
      delete [] lcusty_;
   }
   nlcust_ = 0;
   lcustx_ = lcusty_ = 0;   
   nlcust_ = n+2;
   lcustx_ = new float[nlcust_];
   lcusty_ = new float[nlcust_];
   
   int idx,i;
   idx = 0;
   lcustx_[idx] = 0;
   lcusty_[idx++] = 0;
   for (i=0; i<n; i++)
   {
      lcustx_[idx] = ((float) x[i])/sx;
      lcusty_[idx++] = (float)(1.0-((float) y[i])/sy);
   }   
   lcustx_[idx] = 0;
   lcusty_[idx++] = 0;
   /*bgLog(" hyst low custom x: ");
   for (i=0; i<=n; i++)
      bgLog(" %f", lcustx_[i]);
   bgLog("\n");
   bgLog(" low custom y: ");
   for (i=0; i<=n; i++)
      bgLog(" %f", lcusty_[i]);
   bgLog("\n");*/
}

// 设置低阈值
void BgEdgeDetect::SetCustomLow(double* x, double* y, int n)
{
   if(nlcust_>0)
   {
      delete [] lcustx_;
      delete [] lcusty_;
   }
   nlcust_ = 0;
   lcustx_ = lcusty_ = 0;   
   nlcust_ = n+2;
   lcustx_ = new float[nlcust_];
   lcusty_ = new float[nlcust_];
   
   int idx,i;
   idx = 0;
   lcustx_[idx] = 0;
   lcusty_[idx++] = 0;
   for (i=0; i<n; i++)
   {
      lcustx_[idx] = (float) x[i];
      lcusty_[idx++] = (float) y[i];
   }   
   lcustx_[idx] = 0;
   lcusty_[idx++] = 0;
   /*bgLog(" hyst low custom x: ");
   for (i=0; i<=n; i++)
      bgLog(" %f", lcustx_[i]);
   bgLog("\n");
   bgLog(" low custom y: ");
   for (i=0; i<=n; i++)
      bgLog(" %f", lcusty_[i]);
   bgLog("\n");   */
}
