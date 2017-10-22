// RAList.h: interface for the RAList class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RALIST_H__EC1048FF_2A0F_4AF5_B82E_9428252E01A5__INCLUDED_)
#define AFX_RALIST_H__EC1048FF_2A0F_4AF5_B82E_9428252E01A5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//define Region Adjacency List class prototype
class RAList {

public:

	//============================
	// *** Public Data Members ***
	//============================
	
	////////////RAM Label//////////
	int		label;
	
	////////////RAM Weight/////////
	float	edgeStrength;
	int		edgePixelCount;
	
	////////////RAM Link///////////
	RAList	*next;
	
	//=======================
	// *** Public Methods ***
	//=======================
	
	/*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
	/* Class Constructor and Destructor */
	/*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/
	
	//***Class Constrcutor***
	RAList( void );
	
	//***Class Destructor***
	~RAList( void );
	
	/*/\/\/\/\/\/\/\/\/\/\/\/\*/
	/*  RAM List Manipulation */
	/*\/\/\/\/\/\/\/\/\/\/\/\/*/
	
	//Usage: Insert(entry)
	int Insert(RAList*);		//Insert a region node into the region adjecency list
	
private:
	
	//=============================
	// *** Private Data Members ***
	//=============================
	
	///////current and previous pointer/////
	RAList	*cur, *prev;  // 当前指针，或前一个指针
	
	////////flag///////////
	unsigned char exists;
	
};

#endif
