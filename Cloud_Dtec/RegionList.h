// RegionList.h: interface for the RegionList class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_REGIONLIST_H__69BD941B_8C8C_40EA_A1BB_4547E83DC19B__INCLUDED_)
#define AFX_REGIONLIST_H__69BD941B_8C8C_40EA_A1BB_4547E83DC19B__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//include global type definitions
#include	"tdef.h"

//define region structure
struct REGION {
	int			label;		     // Ψһ��ʶһ�������
	int			boundaryPtCount; // �߽����ص����
	int         regionPtCount;   // ���������ظ���
	int			boundaryPtIdx;   // �߽���������
	int         regionPtIdx;     // ������������
};

//region class prototype...
class RegionList {

public:

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /* Class Constructor and Destructor */
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|               *  Class Constructor  *              |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Constructs a region list object.                 |//
  //|                                                    |//
  //|   Its arguments are:                               |//
  //|                                                    |//
  //|   <* maxRegions *>                                 |//
  //|   The maximum amount of regions that can be class- |//
  //|   ified by the region list.                        |//
  //|                                                    |//
  //|   <* L *>                                          |//
  //|   The length of the input data set being class-    |//
  //|   ified by the region list object.                 |//
  //|                                                    |//
  //|   <* N *>                                          |//
  //|   The dimension of the input data set being class- |//
  //|   ified by the region list object.                 |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|     RegionList(maxRegions, C, N, L)                   |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

	RegionList(int, int, int, int);

	// Class Destructor
	~RegionList( void );

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*  Region List Manipulation  */
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|                *  Add Region  *                    |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Adds a region to the region list.                |//
  //|                                                    |//
  //|   Its arguments are:                               |//
  //|                                                    |//
  //|   <* label *>                                      |//
  //|                                                    |//
  //|   A positive integer used to uniquely identify     |//
  //|   a region.                                        |//
  //|                                                    |//
  //|   <* pointCount *>                                 |//
  //|   A positive integer that specifies the number of  |//
  //|   N-dimensional data points that exist in the re-  |//
  //|   gion being classified.                           |//
  //|                                                    |//
  //|   <* indeces *>                                    |//
  //|   An integer array that specifies the set of ind-  |//
  //|   eces of the data points that are contianed with- |//
  //|   in this region.                                  |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|     AddRegion(label, pointCount, indeces)          |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

	void AddRegion(int, int, int*, int, int*);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|                    *  Reset  *                     |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Resets the region list for re-use (for new       |//
  //|   classification).                                 |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

	void Reset( void );	

  /*/\/\/\/\/\/\/\/\/\/\*/
  /*  Query Region List */
  /*\/\/\/\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|          *  Get Number of Regions  *               |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Returns the number of regions stored by the      |//
  //|   region list.                                     |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

	int	GetNumRegions ( void );

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|                  *  Get Label  *                   |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Returns the label of a specified region.         |//
  //|                                                    |//
  //|   Its arguments are:                               |//
  //|                                                    |//
  //|   <* regionNumber *>                               |//
  //|   The index of the region in the region list       |//
  //|   array.                                           |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|     label = GetLabel(regionNumber)                 |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

	int	GetLabel(int);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|                *  Get Region Count  *              |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Returns number of data points contained by a sp- |//
  //|   ecified region.                                  |//
  //|                                                    |//
  //|   Its arguments are:                               |//
  //|                                                    |//
  //|   <* regionNumber *>                               |//
  //|   The index of the region in the region list       |//
  //|   array.                                           |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|     pointCount = GetRegionCount(regionNumber)      |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

	int GetBourdaryCount(int);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|               *  Get Region Indeces  *             |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Returns a pointer to a set of grid location ind- |//
  //|   eces specifying the data points belonging to a   |//
  //|   specified region.                                |//
  //|                                                    |//
  //|   Its arguments are:                               |//
  //|                                                    |//
  //|   <* regionNumber *>                               |//
  //|   The index of the region in the region list       |//
  //|   array.                                           |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|     indeces = GetBoundaryIndeces(regionNumber)       |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

	int*GetBoundaryIndeces(int regionLable); //// ��ȡregionLable�߽�����

	int GetRegionPtCount(int regionLable); // ��ȡ��ʶΪregionLable����ĵ���
	int*GetRegionPtIndeces(int regionLable); // ��ȡ��ʶΪregionLable���������������
	int GetTotalBoundaryPtCount(void); // ��ȡ�ܵı߽����
	int GetTotalRegionPtCount(void);   // ��ȡ�ܵ��������


private:

  /*/\/\/\/\/\/\/\/\/\/\/\*/
  /*  Class Error Handler */
  /*\/\/\/\/\/\/\/\/\/\/\/*/

	void ErrorHandler(char*, char*, ErrorType);

  //=============================
  // *** Private Data Members ***
  //=============================

	//#####################################
	//### REGION LIST PARTITIONED ARRAY ###
	//#####################################

	REGION		*regionList;			//array of maxRegions regions
	int			minRegion;

	int			maxRegions;				//defines the number maximum number of regions
										//allowed (determined by user during class construction)
	int			numRegions;				//the number of regions currently stored by the
										//region list
	int			freeRegion;				//an index into the regionList pointing to the next
										//available region in the regionList

	//#####################################
	//###         INDEX TABLE           ###
	//#####################################

	int			*boundariesIdxTable;			//an array of indexes that point into an external structure
										//specifying which points belong to a region
	
	int			freeBlockLoc;			//points to the next free block of memory in the indexTable

	//#####################################
	//###     INPUT DATA PARAMETERS     ###
	//#####################################

	//Dimension of data set
	int			N;						//dimension of data set being classified by region list
										//class

	//Length of the data set
	int			C;						//number of points contained by the data set being classified by
										//region list class

	int         L;                      //���������ص���

	int         *regionsIdxTable;   
	int         freeRegionBlockLoc;

};

#endif
