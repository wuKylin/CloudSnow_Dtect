// BgEdgeList.h: interface for the BgEdgeList class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BGEDGELIST_H__8926E172_7715_463C_842D_2EDD592AB980__INCLUDED_)
#define AFX_BGEDGELIST_H__8926E172_7715_463C_842D_2EDD592AB980__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "BgEdge.h"
#include "BgImage.h"

class BgEdgeList

{

public:

	int nEdges_;  // ±ßÔµÊý 
	
	BgEdge* edgelist_; // ±ß±í
	
	BgEdge* crtedge_;  // 
	
	
	
	BgEdgeList();
	
	~BgEdgeList();
	
	void AddEdge(float*, int);
	
	void AddEdge(int*, int nPoints);
	
	void RemoveShortEdges(int);
	
	void SetBinImage(BgImage*);
	
	bool SaveEdgeList(char*);
	
	void SetGradient(float*, float*, float*, int);
	
	void SetNoMark(void);
	
	void GetAllEdgePoints(int*, int*, int*);
	
	
	
};

#endif // !defined(AFX_BGEDGELIST_H__8926E172_7715_463C_842D_2EDD592AB980__INCLUDED_)
