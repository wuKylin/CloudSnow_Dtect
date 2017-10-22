// RegionList.cpp: implementation of the RegionList class.
//
//////////////////////////////////////////////////////////////////////

#include "RegionList.h"
#include "RegionList.h"
#include <stdio.h>
#include <stdlib.h>



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

/*******************************************************/
/*Constructor                                          */
/*******************************************************/
/*Constructor                                          */
/*******************************************************/
/*Pre:                                                 */
/*      - modesPtr is a pointer to an array of modes   */
/*      - maxRegions_ is the maximum number of regions */
/*        that can be defined                          */
/*      - L_ is the number of data points being class- */
/*        ified by the region list class               */
/*      - N is the dimension of the data set being cl- */
/*        assified by the region list class            */
/*Post:                                                */
/*      - a region list object has been properly init- */
/*        ialized.                                     */
/*******************************************************/

RegionList::RegionList(int maxRegions_, int C_, int N_, int L_ )
{

	//Obtain maximum number of regions that can be
	//defined by user
	if((maxRegions = maxRegions_) <= 0)
		ErrorHandler("RegionList", "Maximum number of regions is zero or negative.", FATAL);

	//Obtain dimension of data set being classified by
	//region list class
	if((N = N_) <= 0)
		ErrorHandler("RegionList", "Dimension is zero or negative.", FATAL);

	//Obtain length of input data set...
	if((C = C_) <= 0)
		ErrorHandler("RegionList", "Length of data set is zero or negative.", FATAL);

	if((L = L_) <= 0)
		ErrorHandler("RegionList", "Length of data set is zero or negative.", FATAL);

	//Allocate memory for index table
	if(!(boundariesIdxTable = new int [C]))
		ErrorHandler("RegionList", "Not enough memory.", FATAL);

	if(!(regionsIdxTable = new int [L]))
		ErrorHandler("RegionList", "Not enough memory.", FATAL);

	//Allocate memory for region list array
	if(!(regionList = new REGION [maxRegions]))
		ErrorHandler("RegionList", "Not enough memory.", FATAL);

	//Initialize region list...
	numRegions		= freeRegion = 0;

	//Initialize boundariesIdxTable
	freeBlockLoc	= 0;
	freeRegionBlockLoc = 0;

	//done.
	return;

}

/*******************************************************/
/*Destructor                                           */
/*******************************************************/
/*Destroys region list object.                         */
/*******************************************************/
/*Post:                                                */
/*      - region list object has been properly dest-   */
/*        oyed.                                        */
/*******************************************************/

RegionList::~RegionList( void )
{
	//de-allocate memory...
	delete [] regionList;
	delete [] boundariesIdxTable;

	//done.
	return;
}

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /***  Region List Manipulation  ***/
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Add Region                                           */
/*******************************************************/
/*Adds a region to the region list.                    */
/*******************************************************/
/*Pre:                                                 */
/*      - label is a positive integer used to uniquely */
/*        identify a region                            */
/*      - boundaryPtCount is the number of N-dimensional    */
/*        data points that exist in the region being   */
/*        classified.                                  */
/*      - indeces is a set of indeces specifying the   */
/*        data points contained within this region     */
/*      - boundaryPtCount must be > 0                       */
/*Post:                                                */
/*      - a new region labeled using label and contai- */
/*        ning boundaryPtCount number of points has been    */
/*        added to the region list.                    */
/*******************************************************/

void RegionList::AddRegion(int label, int boundaryPtCount, int *boundaryIndeces ,int regionPtCount, int *regionIndeces )
{

	//make sure that there is enough room for this new region 
	//in the region list array...
	if(numRegions >= maxRegions)
		ErrorHandler("AddRegion", "Not enough memory allocated.", FATAL);

	//make sure that label is positive and point Count > 0...
	if((label < 0)||(boundaryPtCount <= 0)||(regionPtCount<=0))
		ErrorHandler("AddRegion", "Label is negative or number of points in region is invalid.", FATAL);

	//make sure that there is enough memory in the boundariesIdxTable
	//for this region...
	if((freeBlockLoc + boundaryPtCount) > C)
		ErrorHandler("AddRegion", "Adding more points than what is contained in data set.", FATAL);

	if((freeRegionBlockLoc + regionPtCount) > L)
		ErrorHandler("AddRegion", "Adding more points than what is contained in data set.", FATAL);

	//place new region into region list array using
	//freeRegion index
	regionList[freeRegion].label				= label;
	regionList[freeRegion].boundaryPtCount		= boundaryPtCount;
	regionList[freeRegion].boundaryPtIdx		= freeBlockLoc;
	regionList[freeRegion].regionPtCount		= regionPtCount;
	regionList[freeRegion].regionPtIdx			= freeRegionBlockLoc;

	//copy indeces into boundariesIdxTable using freeBlock...
	int i;
	for(i = 0; i < boundaryPtCount; i++)
		boundariesIdxTable[freeBlockLoc+i] = boundaryIndeces[i];

	for(i = 0; i < regionPtCount; i++)
		regionsIdxTable[freeRegionBlockLoc+i] = regionIndeces[i];

	//increment freeBlock to point to the next free
	//block
	freeBlockLoc	+= boundaryPtCount;
	freeRegionBlockLoc	+= regionPtCount;

	//increment freeRegion to point to the next free region
	//also, increment numRegions to indicate that another
	//region has been added to the region list
	freeRegion++;
	numRegions++;

	//done.
	return;

}

/*******************************************************/
/*Reset                                                */
/*******************************************************/
/*Resets the region list.                              */
/*******************************************************/
/*Post:                                                */
/*      - the region list has been reset.              */
/*******************************************************/

void RegionList::Reset( void )
{

	//reset region list
	freeRegion = numRegions = freeBlockLoc = 0;
	freeRegionBlockLoc = 0;

	//done.
	return;

}

  /*/\/\/\/\/\/\/\/\/\/\*/
  /*  Query Region List */
  /*\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Get Number Regions                                   */
/*******************************************************/
/*Returns the number of regions stored by region list. */
/*******************************************************/
/*Post:                                                */
/*      - the number of regions stored by the region   */
/*        list is returned.                            */
/*******************************************************/

int	RegionList::GetNumRegions( void )
{
	// return region count
	return numRegions;
}

/*******************************************************/
/*Get Label                                            */
/*******************************************************/
/*Returns the label of a specified region.             */
/*******************************************************/
/*Pre:                                                 */
/*      - regionNum is an index into the region list   */
/*        array.                                       */
/*Post:                                                */
/*      - the label of the region having region index  */
/*        specified by regionNum has been returned.    */
/*******************************************************/

int	RegionList::GetLabel(int regionNum)
{
	//return the label of a specified region
	return regionList[regionNum].label;
}

/*******************************************************/
/*Get Region Count                                     */
/*******************************************************/
/*Returns the point count of a specified region.       */
/*******************************************************/
/*Pre:                                                 */
/*      - regionNum is an index into the region list   */
/*        array.                                       */
/*Post:                                                */
/*      - the number of points that classify the       */
/*        region whose index is specified by regionNum */
/*        is returned.                                 */
/*******************************************************/
int RegionList::GetRegionPtCount(int regionNum)
{
	return regionList[regionNum].regionPtCount;
}
// 获取周长
int RegionList::GetBourdaryCount(int regionNum)
{
	//return the region count of a specified region
	return regionList[regionNum].boundaryPtCount;
}

/*******************************************************/
/*Get Region Indeces                                   */
/*******************************************************/
/*Returns the point indeces specifying a region.       */
/*******************************************************/
/*Pre:                                                 */
/*      - regionNum is an index into the region list   */
/*        array.                                       */
/*Post:                                                */
/*      - the region indeces specifying the points     */
/*        contained by the region specified by region- */
/*        Num are returned.                            */
/*******************************************************/
int *RegionList::GetRegionPtIndeces(int regionNum)
{
	return &regionsIdxTable[regionList[regionNum].regionPtIdx];
}
int *RegionList::GetBoundaryIndeces(int regionNum)
{
	//return point indeces using regionNum
	return &boundariesIdxTable[regionList[regionNum].boundaryPtIdx];
	// 原来boundariesIdxTable是用来存储每个区域中的所包含的边界像素点的~~
}
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@     PRIVATE METHODS     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

  /*/\/\/\/\/\/\/\/\/\/\/\*/
  /*  Class Error Handler */
  /*\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Error Handler                                        */
/*******************************************************/
/*Class error handler.                                 */
/*******************************************************/
/*Pre:                                                 */
/*      - functName is the name of the function that   */
/*        caused an error                              */
/*      - errmsg is the error message given by the     */
/*        calling function                             */
/*      - status is the error status: FATAL or NON-    */
/*        FATAL                                        */
/*Post:                                                */
/*      - the error message errmsg is flagged on beh-  */
/*        ave of function functName.                   */
/*      - if the error status is FATAL then the program*/
/*        is halted, otherwise execution is continued, */
/*        error recovery is assumed to be handled by   */
/*        the calling function.                        */
/*******************************************************/

void RegionList::ErrorHandler(char *functName, char* errmsg, ErrorType status)
{

	//flag error message on behalf of calling function, error format
	//specified by the error status...
	if(status == NONFATAL)
		fprintf(stderr, "\n%s Error: %s\n", functName, errmsg);
	else
	{
		fprintf(stderr, "\n%s Fatal Error: %s\n\nAborting Program.\n\n", functName, errmsg);
		exit(1);
	}

}


int RegionList::GetTotalBoundaryPtCount(void) // 获取总的边界点数
{
	return C;

}

int RegionList::GetTotalRegionPtCount(void)   // 获取总的区域点数
{
	return L;

}