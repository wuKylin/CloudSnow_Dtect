// BgEdge.h: interface for the BgEdge class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BGEDGE_H__40B461B5_A37D_4672_8277_92AFF7150E21__INCLUDED_)
#define AFX_BGEDGE_H__40B461B5_A37D_4672_8277_92AFF7150E21__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

// 8¡⁄”ÚÕº
// 8-connected neighbour
static const int gNb8[8][2]=
{
   1, 1,
   1,-1,
   1, 0,
   0, 1,
   0,-1,
  -1,-1,
  -1, 0,
  -1, 1
};

// …Ë÷√±ﬂ‘µÕº
class BgEdge
{
public:
   int* edge_;
   double* grad_;
   bool isGradSet_;
   bool isMarkSet_;
   unsigned char* mark_;
   int nPoints_;
   BgEdge* next_;

   BgEdge();
   ~BgEdge();
   void SetPoints(float*, int);
   void SetPoints(int*, int);
   void SetGradient(float*,float*,float*,int);
};

#endif // !defined(AFX_BGEDGE_H__40B461B5_A37D_4672_8277_92AFF7150E21__INCLUDED_)
