// msImageProcessor.cpp: implementation of the msImageProcessor class.
//
//////////////////////////////////////////////////////////////////////


#include "msImageProcessor.h"
#include	<math.h>
#include	<stdio.h>
#include	<assert.h>
#include	<string.h>
#include	<stdlib.h>

//#ifdef _DEBUG
//#undef THIS_FILE
//static char THIS_FILE[]=__FILE__;
//#define new DEBUG_NEW
//#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

/*******************************************************/
/*Class Constructor                                    */
/*******************************************************/
/*Post:                                                */
/*      The msImageProcessor class has been properly   */
/*      initialized.                                   */
/*******************************************************/

msImageProcessor::msImageProcessor( void )
{

	//intialize basin of attraction structure
	//used by the filtering algorithm
	modeTable			= NULL;
	pointList			= NULL;
	pointCount			= 0;

	//initialize region list
	regionList			= NULL;

	//initialize output structures...
	msRawData			= NULL;
	labels				= NULL;
	modes				= NULL;
	modePointCounts		= NULL;
	regionCount			= 0;

	//intialize temporary buffers used for
	//performing connected components
	indexTable			= NULL;
	LUV_data			= NULL;
	// 

	//initialize region adjacency matrix
	raList				= NULL;
	freeRAList			= NULL;
	raPool				= NULL;

	//intialize visit table to having NULL entries
	visitTable			= NULL;

	//initialize epsilon such that transitive closure
	//does not take edge strength into consideration when
	//fusing regions of similar color
	epsilon				= 1.0;

	//initialize class state to indicate that
	//an output data structure has not yet been
	//created...
	class_state.OUTPUT_DEFINED	= false;

    //Changed by Sushil from 1.0 to 0.1, 11/11/2008
    LUV_treshold = 0.1;
}

/*******************************************************/
/*Class Destructor                                     */
/*******************************************************/
/*Post:                                                */
/*      The msImageProcessor class has been properly   */
/*      destroyed.                                     */
/*******************************************************/

msImageProcessor::~msImageProcessor( void )
{

	//de-allocate memory
	if(class_state.OUTPUT_DEFINED)	DestroyOutput();
	if(regionList)					delete regionList;
	regionList = NULL;

	//done.

}

 /*/\/\/\/\/\/\/\/\/\/\/\/\/\*/
 /*  Input Image Declaration */
 /*\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Define Image                                         */
/*******************************************************/
/*Uploads an image into the image segmenter class to   */
/*be segmented.                                        */
/*******************************************************/
/*Pre:                                                 */
/*      - data_ is a one dimensional array of unsigned */
/*        char RGB vectors                             */
/*      - type is the type of the image: COLOR or      */
/*        GREYSCALE                                    */
/*      - height_ and width_ define the dimension of   */
/*        the image                                    */
/*      - if the image is of type GREYSCALE then       */
/*        data containes only one number per pixel     */
/*        location, where a pixel location is defined  */
/*        by the index into the data array             */
/*Post:                                                */
/*      - the image specified has been uploaded into   */
/*        the image segmenter class to be segmented.   */
/*******************************************************/

void msImageProcessor::DefineImage(byte *data_, imageType type, int height_, int width_)
{

	//obtain image dimension from image type
	int dim;
	if(type == COLOR)
		dim	= 3;
	else
		dim = 1;

	//perfor rgb to luv conversion
	int		i;
	float	*luv	= new float [height_*width_*dim];
	if(dim == 1)
	{
		for(i = 0; i < height_*width_; i++)
			luv[i]	= (float)(data_[i]);
	}
	else
	{
		for(i = 0; i < height_*width_; i++)
		{
				RGBtoLUV(&data_[dim*i], &luv[dim*i]);
		}
	}

	//define input defined on a lattice using mean shift base class
	DefineLInput(luv, height_, width_, dim); // 定义输入

	//Define a default kernel if it has not been already
	//defined by user
	if(!h)
	{
		//define default kernel paramerters...
		kernelType	k[2]		= {Uniform, Uniform};
		int			P[2]		= {2, N};
		float		tempH[2]	= {1.0 , 1.0};

		//define default kernel in mean shift base class
		DefineKernel(k, tempH, P, 2);// 定义核函数
	}

	//de-allocate memory
	delete [] luv;

	//done.
	return;

}

void msImageProcessor::DefineBgImage(byte* data_, imageType type, int height_, int width_)
{

	//obtain image dimension from image type
	int dim;
	if(type == COLOR)
		dim	= 3;
	else
		dim = 1;

	//perform texton classification
	int		i;

	float	*luv	= new float [height_*width_*dim];
	
	if(dim == 1)
	{
		for(i = 0; i < height_*width_; i++)
			luv[i]	= (float)(data_[i]);
	}
	else
	{
		for(i = 0; i < height_*width_; i++)
			RGBtoLUV(&data_[dim*i], &luv[dim*i]);
		
	}

	//define input defined on a lattice using mean shift base class
	DefineLInput(luv, height_, width_, dim);


	//Define a default kernel if it has not been already
	//defined by user
	if(!h)
	{
		//define default kernel paramerters...
		kernelType	k[2]		= {Uniform, Uniform};
		int			P[2]		= {2, N};
		float		tempH[2]	= {1.0 , 1.0};

		//define default kernel in mean shift base class
		DefineKernel(k, tempH, P, 2);
	}

	//de-allocate memory
	delete [] luv;

	//done.
	return;

}

 /*/\/\/\/\/\/\/\/\*/
 /*   Weight Map   */
 /*\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Set Weight Map                                       */
/*******************************************************/
/*Populates the weight map with specified edge         */
/*strengths.                                           */
/*******************************************************/
/*Pre:                                                 */
/*      - wm is a floating point array of size         */
/*        (height x width) specifying for each pixel   */
/*        edge strength.                               */
/*      - eps is a threshold used to fuse similar      */
/*        regions during transitive closure.           */
/*Post:                                                */
/*      - wm has been used to populate the weight      */
/*        map.                                         */
/*      - the threshold used during transitive closure */
/*        is taken as eps.                             */
/*******************************************************/

void msImageProcessor::SetWeightMap(float *wm, float eps)
{

	//initlaize confmap using wm
	SetLatticeWeightMap(wm);

	//set threshold value
	if((epsilon = eps) < 0)
		ErrorHandler("msImageProcessor", "SetWeightMap", "Threshold is negative.");

	//done.
	return;

}

/*******************************************************/
/*Remove Weight Map                                    */
/*******************************************************/
/*Removes the weight map.                              */
/*******************************************************/
/*Post:                                                */
/*      - the weight map has been removed.             */
/*      - if a weight map did not exist NO error       */
/*        is flagged.                                  */
/*******************************************************/

void msImageProcessor::RemoveWeightMap( void )
{

	//remove confmap
	RemoveLatticeWeightMap();

	//set threshold value to zero
	epsilon	= 0;

	//done.
	return;

}

 /*/\/\/\/\/\/\/\/\/\*/
 /* Image Filtering  */
 /*\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Filter                                               */
/*******************************************************/
/*Performs mean shift filtering on the specified input */
/*image using a user defined kernel.                   */
/*******************************************************/
/*Pre:                                                 */
/*      - the user defined kernel used to apply mean   */
/*        shift filtering to the defined input image   */
/*        has spatial bandwidth sigmaS and range band- */
/*        width sigmaR                                 */
/*      - speedUpLevel determines whether or not the   */
/*        filtering should be optimized for faster     */
/*        execution: a value of NO_SPEEDUP turns this  */
/*        optimization off and a value SPEEDUP turns   */
/*        this optimization on                         */
/*      - a data set has been defined                  */
/*      - the height and width of the lattice has been */
/*        specified using method DefineLattice()       */
/*Post:                                                */
/*      - mean shift filtering has been applied to the */
/*        input image using a user defined kernel      */
/*      - the filtered image is stored in the private  */
/*        data members of the msImageProcessor class.  */
/*******************************************************/

// 平滑，区域生长
void msImageProcessor::Filter(int sigmaS, float sigmaR, SpeedUpLevel speedUpLevel)
{

	//Check Class consistency...

	//check:
	// (1) if this operation is consistent
	// (2) if kernel was created
	// (3) if data set is defined
	// (4) if the dimension of the kernel agrees with that
	//     of the defined data set
	// if not ... flag an error!
	classConsistencyCheck(N+2, true);
	if(ErrorStatus == EL_ERROR)
		return;

	//If the algorithm has been halted, then exit
	if((ErrorStatus = msSys.Progress((float)(0.0))) == EL_HALT)
	{
		return;
	}
	
	//If the image has just been read then allocate memory
	//for and initialize output data structure used to store
	//image modes and their corresponding regions...
	if(class_state.OUTPUT_DEFINED == false)
	{
		InitializeOutput();

		//check for errors...
		if(ErrorStatus == EL_ERROR)
			return;
	}

	//****************** Allocate Memory ******************

	//Allocate memory for basin of attraction mode structure...
	if((!(modeTable = new unsigned char [L]))||(!(pointList = new int [L])))
	{
		ErrorHandler("msImageProcessor", "Allocate", "Not enough memory.");
		return;
	}

	//start timer
#ifdef PROMPT
	double timer;
	msSys.StartTimer();
#endif

	//*****************************************************

	/*****************************************************根据speedUpLevel进行不同程度的平滑************************************/
	//filter image according to speedup level...
	switch(speedUpLevel)
	{
	//no speedup...
	case NO_SPEEDUP:	
      //NonOptimizedFilter((float)(sigmaS), sigmaR);	break;
      NewNonOptimizedFilter((float)(sigmaS), sigmaR);	
	  break;
	//medium speedup
	case MED_SPEEDUP:	
      //OptimizedFilter1((float)(sigmaS), sigmaR);		break;
      NewOptimizedFilter1((float)(sigmaS), sigmaR);		
	  break;
	//high speedup
	case HIGH_SPEEDUP: 
      //OptimizedFilter2((float)(sigmaS), sigmaR);		break;
      NewOptimizedFilter2((float)(sigmaS), sigmaR);		
	  break;
   // new speedup
	}

	//****************** Deallocate Memory ******************

	//de-allocate memory used by basin of attraction mode structure
	delete [] modeTable;
	delete [] pointList;

	//re-initialize structure
	modeTable	= NULL;
	pointList	= NULL;
	pointCount	= 0;

	//*******************************************************

	//If the algorithm has been halted, then de-allocate the output
	//and exit
	if((ErrorStatus = msSys.Progress((float)(0.8))) == EL_HALT)
	{
		DestroyOutput();
		return;
	}

	//Label image regions, also if segmentation is not to be
	//performed use the resulting classification structure to
	//calculate the image boundaries...

   /*
	//copy msRawData into LUV_data, rounding each component of each
	//LUV value stored by msRawData to the nearest integer
	int	i;
	for(i = 0; i < L*N; i++)
	{
		if(msRawData[i] < 0)
			LUV_data[i] = (int)(msRawData[i] - 0.5);
		else
			LUV_data[i] = (int)(msRawData[i] + 0.5);
	}
   */
   int i;
   for (i=0; i<L*N; i++)
   {
      LUV_data[i] = msRawData[i];
   }


#ifdef PROMPT
	timer	= msSys.ElapsedTime();
	msSys.Prompt("(%6.2f sec)\nConnecting regions         ...", timer);
	msSys.StartTimer();
#endif
	

	/****************************************************区域生长Connect************************************/
	//Perform connecting (label image regions) using LUV_data
	Connect();
	
#ifdef PROMPT
	timer	= msSys.ElapsedTime();
	msSys.Prompt("done. (%6.2f seconds, numRegions = %6d)\n", timer, regionCount);
	msSys.StartTimer();
#endif

	//done.
	return;

}

 /*/\/\/\/\/\/\/\/\/\/\/\*/
 /* Image Region Fusing  */
 /*\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Fuse Regions                                         */
/*******************************************************/
/*Fuses the regions of a filtered image.               */
/*******************************************************/
/*Pre:                                                 */
/*      - the range radius is specified by sigmaR      */
/*      - minRegion is the minimum point density that  */
/*        a region may have in the resulting segment-  */
/*        ed image                                     */
/*      - a data set has been defined                  */
/*      - the height and width of the lattice has been */
/*        specified using method DefineLattice()       */
/*Post:                                                */
/*      - the image regions have been fused.           */
/*      - if an result is stored by this class then    */
/*        this result is used as input to this method. */
/*      - if no result is stored by this class,        */
/*        the input image defined by calling the       */
/*        method DefineImage is used.                  */
/*******************************************************/
// 区域生长，合并，修剪，FuseRegions就是比segment少了一个filter，分割可以看作先filter后fuseRegions
void msImageProcessor::FuseRegions(float sigmaS, int minRegion)
{

	//Check Class consistency...

	//check:
	// (1) if this operation is consistent
	// (2) if kernel was created
	// (3) if data set is defined
	// (4) if the dimension of the kernel agrees with that
	//     of the defined data set
	// if not ... flag an error!
	//********************************************************************这个不懂Check Class consistency...
	classConsistencyCheck(N+2, true);
	if(ErrorStatus == EL_ERROR)
		return;

	//Check to see if the algorithm is to be halted, if so then
	//destroy output and exit
	if((ErrorStatus = msSys.Progress((float)(0.8))) == EL_HALT)
	{
		if(class_state.OUTPUT_DEFINED)	DestroyOutput();
		return;
	}

	//obtain sigmaS (make sure it is not zero or negative, if not
	//flag an error)
	if((h[1] = sigmaS) <= 0)
	{
		ErrorHandler("msImageProcessor", "FuseRegions", "The feature radius must be greater than or equal to zero.");
		return;
	}

	//if output has not yet been generated then classify the input
	//image regions to be fused...
	if(!(class_state.OUTPUT_DEFINED))
	{

		//Initialize output data structure used to store
		//image modes and their corresponding regions...
		InitializeOutput();
		
		//check for errors...
		if(ErrorStatus == EL_ERROR)
			return;

		//copy data into LUV_data used to classify
		//image regions
      /*
		int i;
		for(i = 0; i < L*N; i++)
		{
			if(data[i] < 0)
				LUV_data[i] = (int)(data[i] - 0.5);
			else
				LUV_data[i] = (int)(data[i] + 0.5);
		}
      */
      int i;
      for (i=0; i<L*N; i++)
      {
         LUV_data[i] = data[i];
      }
		
#ifdef PROMPT
		msSys.Prompt("Connecting regions         ...");
		msSys.StartTimer();
#endif

		//Perform connecting (label image regions) using LUV_data
	/******************************************************区域生长Connect*************************************************************************/
		Connect();
		
		//check for errors
		if(ErrorStatus == EL_ERROR)
			return;
		
#ifdef PROMPT
		double timer	= msSys.ElapsedTime();
		msSys.Prompt("done. (%6.2f seconds, numRegions = %6d)\n", timer, regionCount);
#endif
		
	}

	//Check to see if the algorithm is to be halted, if so then
	//destroy output and exit
	if((ErrorStatus = msSys.Progress((float)(0.85))) == EL_HALT)
	{
		DestroyOutput();
		return;
	}

#ifdef PROMPT
	msSys.Prompt("Applying transitive closure...");
	msSys.StartTimer();
#endif

	//allocate memory visit table
	visitTable = new unsigned char [L];

	//Apply transitive closure iteratively to the regions classified
	//by the RAM updating labels and modes until the color of each neighboring
	//region is within sqrt(rR2) of one another.
	rR2 = (float)(h[1]*h[1]*0.25);
	/**********************************************************重要的传递闭包 *******************************************/
	TransitiveClosure();
	int oldRC = regionCount;
	int deltaRC, counter = 0;
	do {
		TransitiveClosure();
		deltaRC = oldRC-regionCount;
		oldRC = regionCount;
		counter++;
	} while ((deltaRC <= 0)&&(counter < 10));

	//de-allocate memory for visit table
	delete [] visitTable;
	visitTable	= NULL;

	//Check to see if the algorithm is to be halted, if so then
	//destroy output and region adjacency matrix and exit
	if((ErrorStatus = msSys.Progress((float)(1.0))) == EL_HALT)
	{
		DestroyRAM();
		DestroyOutput();
		return;
	}

#ifdef PROMPT
	double timer	= msSys.ElapsedTime();
	msSys.Prompt("done. (%6.2f seconds, numRegions = %6d)\nPruning spurious regions   ...", timer, regionCount);
	msSys.StartTimer();
#endif
	/*******************************************合并小的多边形prune*****************************************************/
	//Prune spurious regions (regions whose area is under
	//minRegion) using RAM
	Prune(minRegion);

#ifdef PROMPT
	timer	= msSys.ElapsedTime();
	msSys.Prompt("done. (%6.2f seconds, numRegions = %6d)\n", timer, regionCount);
	msSys.StartTimer();
#endif

	//Check to see if the algorithm is to be halted, if so then
	//destroy output and region adjacency matrix and exit
	if((ErrorStatus = msSys.Progress((float)(1.0))) == EL_HALT)
	{
		DestroyRAM();
		DestroyOutput();
		return;
	}

	//de-allocate memory for region adjacency matrix
	DestroyRAM();

	//output to msRawData
	int i, j, label;
	for(i = 0; i < L; i++)
	{
		label	= labels[i];
		for(j = 0; j < N; j++)
			{
				msRawData[N*i+j] = modes[N*label+j];
			}
	}

	//done.
	return;

}

  /*/\/\/\/\/\/\/\/\/\/\*/
  /* Image Segmentation */
  /*\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Segment                                              */
/*******************************************************/
/*Segments the defined image.                          */
/*******************************************************/
/*Pre:                                                 */
/*      - sigmaS and sigmaR are the spatial and range  */
/*        radii of the search window respectively      */
/*      - minRegion is the minimum point density that  */
/*        a region may have in the resulting segment-  */
/*        ed image                                     */
/*      - speedUpLevel determines whether or not the   */
/*        filtering should be optimized for faster     */
/*        execution: a value of NO_SPEEDUP turns this  */
/*        optimization off and a value SPEEDUP turns   */
/*        this optimization on                         */
/*Post:                                                */
/*      - the defined image is segmented and the       */
/*        resulting segmented image is stored in the   */
/*        private data members of the image segmenter  */
/*        class.                                       */
/*      - any regions whose point densities are less   */
/*        than or equal to minRegion have been pruned  */
/*        from the segmented image.                    */
/*******************************************************/

// filter,TransitiveClosure,prune
void msImageProcessor::Segment(int sigmaS, float sigmaR, int minRegion, SpeedUpLevel speedUpLevel)
{

	//make sure kernel is properly defined...
	if((!h)||(kp < 2))
	{
		ErrorHandler("msImageProcessor", "Segment", "Kernel corrupt or undefined.");
		return;
	}

	/************************************重要的图像过滤平滑，有connect()*********************************************/
	//Apply mean shift to data set using sigmaS and sigmaR...
	Filter(sigmaS, sigmaR, speedUpLevel);

	//check for errors
	if(ErrorStatus == EL_ERROR)
		return;

	//check to see if the system has been halted, if so exit
	if(ErrorStatus == EL_HALT)
		return;

	//Check to see if the algorithm is to be halted, if so then
	//destroy output and exit
	if((ErrorStatus = msSys.Progress((float)(0.85))) == EL_HALT)
	{
		DestroyOutput();
		return;
	}

#ifdef PROMPT
	msSys.Prompt("Applying transitive closure...");
	msSys.StartTimer();
#endif

	//allocate memory visit table
	visitTable = new unsigned char [L];

	//Apply transitive closure iteratively to the regions classified
	//by the RAM updating labels and modes until the color of each neighboring
	//region is within sqrt(rR2) of one another.
	rR2 = (float)(h[1]*h[1]*0.25);
	/************************************重要的传递闭包！有BuildRAM和computeEdgeStrengths*********************************************/
	
	TransitiveClosure();
	int oldRC = regionCount;
	int deltaRC, counter = 0;
	
	
	do {
		TransitiveClosure();
		deltaRC = oldRC-regionCount;
		oldRC = regionCount;
		counter++;
	} while ((deltaRC <= 0)&&(counter < 10));

	//de-allocate memory for visit table
	delete [] visitTable;
	visitTable	= NULL;

	//Check to see if the algorithm is to be halted, if so then
	//destroy output and regions adjacency matrix and exit
	if((ErrorStatus = msSys.Progress((float)(0.95))) == EL_HALT)
	{
		DestroyRAM();
		DestroyOutput();
		return;
	}

#ifdef PROMPT
	double timer	= msSys.ElapsedTime();
	msSys.Prompt("done. (%6.2f seconds, numRegions = %6d).\nPruning spurious regions\t... ", timer, regionCount);
	msSys.StartTimer();
#endif

	//Prune spurious regions (regions whose area is under
	//minRegion) using RAM
	/************************************重要的多边形修剪！！,有BuildRAM*********************************************/
	Prune(minRegion);

#ifdef PROMPT
	timer	= msSys.ElapsedTime();
	msSys.Prompt("done. (%6.2f seconds, numRegions = %6d)\nPruning spurious regions    ...", timer, regionCount);
	msSys.StartTimer();
#endif

	//Check to see if the algorithm is to be halted, if so then
	//destroy output and regions adjacency matrix and exit
	if((ErrorStatus = msSys.Progress(1.0)) == EL_HALT)
	{
		DestroyRAM();
		DestroyOutput();
		return;
	}

	//de-allocate memory for region adjacency matrix
	DestroyRAM();

	//output to msRawData
	int j, i, label;
	for(i = 0; i < L; i++)
	{
		label	= labels[i];
		for(j = 0; j < N; j++)
			{
				msRawData[N*i+j] = modes[N*label+j];
			}
	}

	//done.
	return;

}

  /*/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*  Data Space Conversion */
  /*\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*RGB To LUV                                           */
/*******************************************************/
/*Converts an RGB vector to LUV.                       */
/*                                                     */
/*See:                                                 */
/*   G. Wyszecki and W.S. Stiles: Color Science:       */
/*   Concepts and Methods, Quantitative Data and       */
/*   Formulae, Wiley, New York, 1982.                  */
/*******************************************************/
/*Pre:                                                 */
/*      - rgbVal is an unsigned char array containing  */
/*        the RGB vector                               */
/*      - luvVal is a floating point array containing  */
/*        the resulting LUV vector                     */
/*Post:                                                */
/*      - rgbVal has been converted to LUV and the     */
/*        result has been stored in luvVal.            */
/*******************************************************/

void msImageProcessor::RGBtoLUV(byte *rgbVal, float *luvVal)
{

	//delcare variables
	double	x, y, z, L0, u_prime, v_prime, constant;

	//convert RGB to XYZ...
	x		= XYZ[0][0]*rgbVal[0] + XYZ[0][1]*rgbVal[1] + XYZ[0][2]*rgbVal[2];
	y		= XYZ[1][0]*rgbVal[0] + XYZ[1][1]*rgbVal[1] + XYZ[1][2]*rgbVal[2];
	z		= XYZ[2][0]*rgbVal[0] + XYZ[2][1]*rgbVal[1] + XYZ[2][2]*rgbVal[2];

	//convert XYZ to LUV...

	//compute L*
	L0		= y / (255.0 * Yn);
	if(L0 > Lt)
		luvVal[0]	= (float)(116.0 * (pow(L0, 1.0/3.0)) - 16.0);
	else
		luvVal[0]	= (float)(903.3 * L0);

	//compute u_prime and v_prime
	constant	= x + 15 * y + 3 * z;
	if(constant != 0)
	{
		u_prime	= (4 * x) / constant;
		v_prime = (9 * y) / constant;
	}
	else
	{
		u_prime	= 4.0;
		v_prime	= 9.0/15.0;
	}

	//compute u* and v*
    luvVal[1] = (float) (13 * luvVal[0] * (u_prime - Un_prime));
    luvVal[2] = (float) (13 * luvVal[0] * (v_prime - Vn_prime));

	//done.
	return;

}

/*******************************************************/
/*LUV To RGB                                           */
/*******************************************************/
/*Converts an LUV vector to RGB.                       */
/*******************************************************/
/*Pre:                                                 */
/*      - luvVal is a floating point array containing  */
/*        the LUV vector                               */
/*      - rgbVal is an unsigned char array containing  */
/*        the resulting RGB vector                     */
/*Post:                                                */
/*      - luvVal has been converted to RGB and the     */
/*        result has been stored in rgbVal.            */
/*******************************************************/

//define inline rounding function...
inline int my_round(double in_x)
{
	if (in_x < 0)
		return (int)(in_x - 0.5);
	else
		return (int)(in_x + 0.5);
}

void msImageProcessor::LUVtoRGB(float *luvVal, byte *rgbVal)
{

	//declare variables...
	int		r, g, b;
	double	x, y, z, u_prime, v_prime;

	//perform conversion
	if(luvVal[0] < 0.1)
		r = g = b = 0;
	else
	{
		//convert luv to xyz...
		if(luvVal[0] < 8.0)
			y	= Yn * luvVal[0] / 903.3;
		else
		{
			y	= (luvVal[0] + 16.0) / 116.0;
			y  *= Yn * y * y;
		}

		u_prime	= luvVal[1] / (13 * luvVal[0]) + Un_prime;
		v_prime	= luvVal[2] / (13 * luvVal[0]) + Vn_prime;

		x		= 9 * u_prime * y / (4 * v_prime);
		z		= (12 - 3 * u_prime - 20 * v_prime) * y / (4 * v_prime);

		//convert xyz to rgb...
		//[r, g, b] = RGB*[x, y, z]*255.0
		r		= my_round((RGB[0][0]*x + RGB[0][1]*y + RGB[0][2]*z)*255.0);
		g		= my_round((RGB[1][0]*x + RGB[1][1]*y + RGB[1][2]*z)*255.0);
		b		= my_round((RGB[2][0]*x + RGB[2][1]*y + RGB[2][2]*z)*255.0);

		//check bounds...
		if(r < 0)	r = 0; if(r > 255)	r = 255;
		if(g < 0)	g = 0; if(g > 255)	g = 255;
		if(b < 0)	b = 0; if(b > 255)	b = 255;

	}

	//assign rgb values to rgb vector rgbVal
	rgbVal[0]	= r;
	rgbVal[1]	= g;
	rgbVal[2]	= b;

	//done.
	return;

}

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*  Filtered and Segmented Image Output */
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Get Raw Data                                         */
/*******************************************************/
/*The output image data is returned.                   */
/*******************************************************/
/*Pre:                                                 */
/*      - outputImageData is a pre-allocated floating  */
/*        point array used to store the filtered or    */
/*        segmented image pixels.                      */
/*Post:                                                */
/*      - the filtered or segmented image data is      */
/*        stored by outputImageData.                   */
/*******************************************************/

void msImageProcessor::GetRawData(float *outputImageData)
{
	//make sure that outputImageData is not NULL
	if(!outputImageData)
	{
		ErrorHandler("msImageProcessor", "GetRawData", "Output image data buffer is NULL.");
		return;
	}

	//copy msRawData to outputImageData
	int i;
	for(i = 0; i < L*N; i++)
		outputImageData[i] = msRawData[i];

	//done.
	return;
}

/*******************************************************/
/*Get Results                                          */
/*******************************************************/
/*The output image is returned.                        */
/*******************************************************/
/*Pre:                                                 */
/*      - outputImage is a pre-allocated unsinged char */
/*        array used to store the filtered or segment- */
/*        ed image pixels                              */
/*Post:                                                */
/*      - the filtered or segmented image is stored by */
/*        outputImage.                                 */
/*******************************************************/

void msImageProcessor::GetResults(byte *outputImage)
{

	//make sure that outpuImage is not NULL
	if(!outputImage)
	{
		ErrorHandler("msImageProcessor", "GetResults", "Output image buffer is NULL.");
		return;
	}

	//if the image type is GREYSCALE simply
	//copy it over to the segmentedImage
	if(N == 1)
	{
		//copy over msRawData to segmentedImage checking
		//bounds
		int i, pxValue;
		for(i = 0; i < L; i++)
		{

			//get value
			pxValue = (int)(msRawData[i]+0.5);
			
			//store into segmented image checking bounds...
			if(pxValue < 0)
				outputImage[i] = (byte)(0);
			else if(pxValue > 255)
				outputImage[i] = (byte)(255);
			else
				outputImage[i] = (byte)(pxValue);

		}

	}
	else if (N == 3)
	{
		
		//otherwise convert msRawData from LUV to RGB
		//storing the result in segmentedImage
		int i;
		for(i = 0; i < L; i++)
			LUVtoRGB(&msRawData[N*i], &outputImage[N*i]);// 

	}
	else
		//Unknown image type: should use MeanShift::GetRawData()...
		ErrorHandler("msImageProcessor", "GetResults", "Unknown image type. Try using MeanShift::GetRawData().");

	//done.
	return;

}

/*******************************************************/
/*Get Boundaries                                       */
/*******************************************************/
/*A region list containing the boundary locations for  */
/*each region is returned.                             */
/*******************************************************/
/*Post:                                                */
/*      - a region list object containing the boundary */
/*        locations for each region is constructed     */
/*      - the region list is returned                  */
/*      - NULL is returned if the image has not been   */
/*        filtered or segmented                        */
/*******************************************************/

RegionList *msImageProcessor::GetBoundariesAndRegions( void )
{

	//define bounds using label information
	if(class_state.OUTPUT_DEFINED)
		DefineRegions();

	//return region list structure
	return regionList;

}
 
/*******************************************************/
/*Get Regions                                          */
/*******************************************************/
/*Returns the regions of the processed image.          */
/*******************************************************/
/*Pre:                                                 */
/*      - labels_out is an integer array of size       */
/*        height*width that stores for each pixel a    */
/*        label relating that pixel to a corresponding */
/*        region in the image                          */
/*      - modes_out is floating point array of size    */
/*        regionCount*N storing the feature component  */
/*        of each region, and indexed by region label  */
/*      - modePointCounts is an integer array of size  */
/*        regionCount, indexed by region label, that   */
/*        stores the area of each region in pixels.    */
/*Post:                                                */
/*      If an input image was defined and processed,   */
/*      - memory has been allocated for labels_out,    */
/*        modes_out and MPC_out.                       */
/*      - labels_out, modes_out, and MPC_out have been */
/*        populated.                                   */
/*      - the number of regions contained by the segm- */
/*        ented image has been returned.               */
/*      If the image has not been defined or processed */
/*      or if there is in-sufficient memory,           */
/*      - no memory has been allocated for labels_out, */
/*        modes_out, and MPC_out.                      */
/*      - -1 is returned for regionCount.              */
/*******************************************************/

int msImageProcessor::GetRegions(int **labels_out, float **modes_out, int **MPC_out)
{
	//check to see if output has been defined for the given input image...
	if(class_state.OUTPUT_DEFINED == false)
		return -1;

	//allocate memory for labels_out, modes_out and MPC_out based
	//on output storage structure
	int		*labels_	= *labels_out, *MPC_out_ = *MPC_out;
	float	*modes_		= *modes_out; 
	if(!(labels_ = new int [L]))
	{
		ErrorHandler("msImageProcessor", "GetRegions", "Not enough memory.");
		return -1;
	}
	if(!(modes_	= new float [regionCount*N]))
	{
		ErrorHandler("msImageProcessor", "GetRegions", "Not enough memory.");
		return -1;
	}
	if(!(MPC_out_ = new int [regionCount]))
	{
		ErrorHandler("msImageProcessor", "GetRegions", "Not enough memory.");
		return -1;
	}

	//populate labels_out with image labels
	int	i;
	for(i = 0; i < L; i++)
		labels_[i] = labels[i];

	//populate modes_out and MPC_out with the color and point
	//count of each region
	for(i = 0; i < regionCount*N; i++)
		modes_[i]	= modes[i];      //存储的是所用数据的模态数据，用label标识一个多边形，同时存储其模态数据，
	//存储每个多边形的模态特征数据
	for(i = 0; i < regionCount; i++)
		MPC_out_[i]	= modePointCounts[i]; //存储每种模式中所包含的点数，一个多边形一个模态，存储多边形的面积

	//done. Return the number of regions resulting from filtering or segmentation.
	return regionCount;
}

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@     PRIVATE METHODS     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

	/*/\/\/\/\/\/\/\/\/\*/
	/*  Image Filtering */
	/*\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Non Optimized Filter                                 */
/*******************************************************/
/*Performs mean shift filtering on the specified input */
/*image using a user defined kernel.                   */
/*******************************************************/
/*Pre:                                                 */
/*      - the user defined kernel used to apply mean   */
/*        shift filtering to the defined input image   */
/*        has spatial bandwidth sigmaS and range band- */
/*        width sigmaR                                 */
/*      - a data set has been defined                  */
/*      - the height and width of the lattice has been */
/*        specified using method DefineLattice()       */
/*Post:                                                */
/*      - mean shift filtering has been applied to the */
/*        input image using a user defined kernel      */
/*      - the filtered image is stored in the private  */
/*        data members of the msImageProcessor class.  */
/*******************************************************/

void msImageProcessor::NonOptimizedFilter(float sigmaS, float sigmaR)
{

	// Declare Variables
	int   iterationCount, i, j;
	double mvAbs;
	
	//make sure that a lattice height and width have
	//been defined...
	if(!height)
	{
		ErrorHandler("msImageProcessor", "LFilter", "Lattice height and width are undefined.");
		return;
	}

	//re-assign bandwidths to sigmaS and sigmaR
	if(((h[0] = sigmaS) <= 0)||((h[1] = sigmaR) <= 0))
	{
		ErrorHandler("msImageProcessor", "Segment", "sigmaS and/or sigmaR is zero or negative.");
		return;
	}
	
	//define input data dimension with lattice
	int lN	= N + 2;
	
	// Traverse each data point applying mean shift
	// to each data point
	
	// Allcocate memory for yk
	double	*yk		= new double [lN];
	
	// Allocate memory for Mh
	double	*Mh		= new double [lN];
	
	// proceed ...
#ifdef PROMPT
	msSys.Prompt("done.\nApplying mean shift (Using Lattice)... ");
#ifdef SHOW_PROGRESS
	msSys.Prompt("\n 0%%");
#endif
#endif

	for(i = 0; i < L; i++)
	{

		// Assign window center (window centers are
		// initialized by createLattice to be the point
		// data[i])
		yk[0] = i%width;
		yk[1] = i/width;
		for(j = 0; j < N; j++)
			yk[j+2] = data[N*i+j];
		
		// Calculate the mean shift vector using the lattice
		LatticeMSVector(Mh, yk);
		
		// Calculate its magnitude squared
		mvAbs = 0;
		for(j = 0; j < lN; j++)
			mvAbs += Mh[j]*Mh[j];
		
		// Keep shifting window center until the magnitude squared of the
		// mean shift vector calculated at the window center location is
		// under a specified threshold (Epsilon)
		
		// NOTE: iteration count is for speed up purposes only - it
		//       does not have any theoretical importance
		iterationCount = 1;
		while((mvAbs >= EPSILON)&&(iterationCount < LIMIT))
		{
			
			// Shift window location
			for(j = 0; j < lN; j++)
				yk[j] += Mh[j];
			
			// Calculate the mean shift vector at the new
			// window location using lattice
			LatticeMSVector(Mh, yk);
			
			// Calculate its magnitude squared
			mvAbs = 0;
			for(j = 0; j < lN; j++)
				mvAbs += Mh[j]*Mh[j];

			// Increment interation count
			iterationCount++;
			
		}

		// Shift window location
		for(j = 0; j < lN; j++)
			yk[j] += Mh[j];
		
		//store result into msRawData...
		for(j = 0; j < N; j++)
			msRawData[N*i+j] = (float)(yk[j+2]);

		// Prompt user on progress
#ifdef SHOW_PROGRESS
		percent_complete = (float)(i/(float)(L))*100;
		msSys.Prompt("\r%2d%%", (int)(percent_complete + 0.5));
#endif
	
		// Check to see if the algorithm has been halted
		if((i%PROGRESS_RATE == 0)&&((ErrorStatus = msSys.Progress((float)(i/(float)(L))*(float)(0.8)))) == EL_HALT)
			break;
	}
	
	// Prompt user that filtering is completed
#ifdef PROMPT
#ifdef SHOW_PROGRESS
	msSys.Prompt("\r");
#endif
	msSys.Prompt("done.");
#endif
	
	// de-allocate memory
	delete [] yk;
	delete [] Mh;

	// done.
	return;

}



/*******************************************************/
/*Optimized Filter 1                                   */
/*******************************************************/
/*Performs mean shift filtering on the specified input */
/*image using a user defined kernel. Previous mode     */
/*information is used to avoid re-applying mean shift  */
/*on certain data points to improve performance.       */
/*******************************************************/
/*Pre:                                                 */
/*      - the user defined kernel used to apply mean   */
/*        shift filtering to the defined input image   */
/*        has spatial bandwidth sigmaS and range band- */
/*        width sigmaR                                 */
/*      - a data set has been defined                  */
/*      - the height and width of the lattice has been */
/*        specified using method DefineLattice()       */
/*Post:                                                */
/*      - mean shift filtering has been applied to the */
/*        input image using a user defined kernel      */
/*      - the filtered image is stored in the private  */
/*        data members of the msImageProcessor class.  */
/*******************************************************/

void msImageProcessor::OptimizedFilter1(float sigmaS, float sigmaR)
{

	// Declare Variables
	int		iterationCount, i, j, k, s, p, modeCandidateX, modeCandidateY, modeCandidate_i;
	float	*modeCandidatePoint;
	double	mvAbs, diff, el;
	
	//make sure that a lattice height and width have
	//been defined...
	if(!height)
	{
		ErrorHandler("msImageProcessor", "LFilter", "Lattice height and width are undefined.");
		return;
	}

	//re-assign bandwidths to sigmaS and sigmaR
	if(((h[0] = sigmaS) <= 0)||((h[1] = sigmaR) <= 0))
	{
		ErrorHandler("msImageProcessor", "Segment", "sigmaS and/or sigmaR is zero or negative.");
		return;
	}
	
	//define input data dimension with lattice
	int lN	= N + 2;
	
	// Traverse each data point applying mean shift
	// to each data point
	
	// Allcocate memory for yk
	double	*yk		= new double [lN];
	
	// Allocate memory for Mh
	double	*Mh		= new double [lN];
	
	// Initialize mode table used for basin of attraction
	memset(modeTable, 0, width*height);

	// Allocate memory mode candidate data point...
	//floating point version
	modeCandidatePoint	= new float	[N];
	
	// proceed ...
#ifdef PROMPT
	msSys.Prompt("done.\nApplying mean shift (Using Lattice) ... ");
#ifdef SHOW_PROGRESS
	msSys.Prompt("\n 0%%");
#endif
#endif


	for(i = 0; i < L; i++)
	{
		// if a mode was already assigned to this data point
		// then skip this point, otherwise proceed to
		// find its mode by applying mean shift...
		if (modeTable[i] == 1)
			continue;

		// initialize point list...
		pointCount = 0;

		// Assign window center (window centers are
		// initialized by createLattice to be the point
		// data[i])
		yk[0] = i%width;
		yk[1] = i/width;
		for(j = 0; j < N; j++)
			yk[j+2] = data[N*i+j];
		
		// Calculate the mean shift vector using the lattice
		LatticeMSVector(Mh, yk);
		
		// Calculate its magnitude squared
		mvAbs = 0;
		for(j = 0; j < lN; j++)
			mvAbs += Mh[j]*Mh[j];
		
		// Keep shifting window center until the magnitude squared of the
		// mean shift vector calculated at the window center location is
		// under a specified threshold (Epsilon)
		
		// NOTE: iteration count is for speed up purposes only - it
		//       does not have any theoretical importance
		iterationCount = 1;
		while((mvAbs >= EPSILON)&&(iterationCount < LIMIT))
		{
			
			// Shift window location
			for(j = 0; j < lN; j++)
				yk[j] += Mh[j];
			
			// check to see if the current mode location is in the
			// basin of attraction...

			// calculate the location of yk on the lattice
			modeCandidateX	= (int) (yk[0]+0.5);
			modeCandidateY	= (int) (yk[1]+0.5);
			modeCandidate_i	= modeCandidateY*width + modeCandidateX;

			// if mvAbs != 0 (yk did indeed move) then check
			// location basin_i in the mode table to see if
			// this data point either:
			
			// (1) has not been associated with a mode yet
			//     (modeTable[basin_i] = 0), so associate
			//     it with this one
			//
			// (2) it has been associated with a mode other
			//     than the one that this data point is converging
			//     to (modeTable[basin_i] = 1), so assign to
			//     this data point the same mode as that of basin_i

			if ((modeTable[modeCandidate_i] != 2) && (modeCandidate_i != i))
			{
				// obtain the data point at basin_i to
				// see if it is within h*TC_DIST_FACTOR of
				// of yk
				for (j = 0; j < N; j++)
					modeCandidatePoint[j] = data[N*modeCandidate_i + j];
				
				// check basin on non-spatial data spaces only
				k = 1;
				s = 0;
				diff = 0;
				while ((diff < TC_DIST_FACTOR) && (k<kp))
				{
					diff = 0;
					for (p=0; p<P[k]; p++)
					{
						el = (modeCandidatePoint[p+s]-yk[p+s+2])/h[k];
						diff += el*el;
					}
					s+=P[k];
					k++;
				}

				// if the data point at basin_i is within
				// a distance of h*TC_DIST_FACTOR of yk
				// then depending on modeTable[basin_i] perform
				// either (1) or (2)
				if (diff < TC_DIST_FACTOR)
				{
					// if the data point at basin_i has not
					// been associated to a mode then associate
					// it with the mode that this one will converge
					// to
					if (modeTable[modeCandidate_i] == 0)
					{
						// no mode associated yet so associate
						// it with this one...
						pointList[pointCount++]		= modeCandidate_i;
						modeTable[modeCandidate_i]	= 2;

					} else
					{

						// the mode has already been associated with
						// another mode, thererfore associate this one
						// mode and the modes in the point list with
						// the mode associated with data[basin_i]...

						// store the mode info into yk using msRawData...
						for (j = 0; j < N; j++)
							yk[j+2] = msRawData[modeCandidate_i*N+j];

						// update mode table for this data point
						// indicating that a mode has been associated
						// with it
						modeTable[i] = 1;

                  // indicate that a mode has been associated
						// to this data point (data[i])
						mvAbs = -1;

						// stop mean shift calculation...
						break;
					}
				}
			}
			
			// Calculate the mean shift vector at the new
			// window location using lattice
			LatticeMSVector(Mh, yk);
			
			// Calculate its magnitude squared
			mvAbs = 0;
			for(j = 0; j < lN; j++)
				mvAbs += Mh[j]*Mh[j];

			// Increment iteration count
			iterationCount++;
			
		}

		// if a mode was not associated with this data point
		// yet associate it with yk...
		if (mvAbs >= 0)
		{
			// Shift window location
			for(j = 0; j < lN; j++)
				yk[j] += Mh[j];
			
			// update mode table for this data point
			// indicating that a mode has been associated
			// with it
			modeTable[i] = 1;
		}
		
		// associate the data point indexed by
		// the point list with the mode stored
		// by yk
		for (j = 0; j < pointCount; j++)
		{
			// obtain the point location from the
			// point list
			modeCandidate_i = pointList[j];

			// update the mode table for this point
			modeTable[modeCandidate_i] = 1;

			//store result into msRawData...
			for(k = 0; k < N; k++)
				msRawData[N*modeCandidate_i+k] = (float)(yk[k+2]);
		}


		//store result into msRawData...
		for(j = 0; j < N; j++)
			msRawData[N*i+j] = (float)(yk[j+2]);

		// Prompt user on progress
#ifdef SHOW_PROGRESS
		percent_complete = (float)(i/(float)(L))*100;
		msSys.Prompt("\r%2d%%", (int)(percent_complete + 0.5));
#endif
	
		// Check to see if the algorithm has been halted
		if((i%PROGRESS_RATE == 0)&&((ErrorStatus = msSys.Progress((float)(i/(float)(L))*(float)(0.8)))) == EL_HALT)
			break;		
	}
	
	// Prompt user that filtering is completed
#ifdef PROMPT
#ifdef SHOW_PROGRESS
	msSys.Prompt("\r");
#endif
	msSys.Prompt("done.");
#endif
	
	// de-allocate memory
	delete [] modeCandidatePoint;
	delete [] yk;
	delete [] Mh;
	
	// done.
	return;

}

/*******************************************************/
/*Optimized Filter 2                                   */
/*******************************************************/
/*Performs mean shift filtering on the specified input */
/*image using a user defined kernel. Previous mode     */
/*information is used to avoid re-applying mean shift  */
/*on certain data points to improve performance. To    */
/*further improve perfmance (during segmentation) poi- */
/*nts within h of a window center during the window    */
/*center's traversal to a mode are associated with the */
/*mode that the window converges to.                   */
/*******************************************************/
/*Pre:                                                 */
/*      - the user defined kernel used to apply mean   */
/*        shift filtering to the defined input image   */
/*        has spatial bandwidth sigmaS and range band- */
/*        width sigmaR                                 */
/*      - a data set has been defined                  */
/*      - the height and width of the lattice has been */
/*        specified using method DefineLattice()       */
/*Post:                                                */
/*      - mean shift filtering has been applied to the */
/*        input image using a user defined kernel      */
/*      - the filtered image is stored in the private  */
/*        data members of the msImageProcessor class.  */
/*******************************************************/

void msImageProcessor::OptimizedFilter2(float sigmaS, float sigmaR)
{

	//if confidence map is null set it to zero
	if(!weightMap)
	{
		weightMap = new float [L];
		int i;
		for(i = 0; i < L; i++)
			weightMap[i] = 0;
	}

	// Declare Variables
	int		iterationCount, i, j, k, s, p, modeCandidateX, modeCandidateY, modeCandidate_i;
	float	*modeCandidatePoint;
	double	mvAbs, diff, el;
	
	//make sure that a lattice height and width have
	//been defined...
	if(!height)
	{
		ErrorHandler("msImageProcessor", "LFilter", "Lattice height and width are undefined.");
		return;
	}

	//re-assign bandwidths to sigmaS and sigmaR
	if(((h[0] = sigmaS) <= 0)||((h[1] = sigmaR) <= 0))
	{
		ErrorHandler("msImageProcessor", "Segment", "sigmaS and/or sigmaR is zero or negative.");
		return;
	}
	
	//define input data dimension with lattice
	int lN	= N + 2;
	
	// Traverse each data point applying mean shift
	// to each data point
	
	// Allcocate memory for yk
	double	*yk		= new double [lN];
	
	// Allocate memory for Mh
	double	*Mh		= new double [lN];
	
	// Initialize mode table used for basin of attraction
	memset(modeTable, 0, width*height);

	// Allocate memory mode candidate data point...
	//floating point version
	modeCandidatePoint	= new float	[N];
	
	// proceed ...
#ifdef PROMPT
	msSys.Prompt("done.\nApplying mean shift (Using Lattice)... ");
#ifdef SHOW_PROGRESS
	msSys.Prompt("\n 0%%");
#endif
#endif

	for(i = 0; i < L; i++)
	{
		// if a mode was already assigned to this data point
		// then skip this point, otherwise proceed to
		// find its mode by applying mean shift...
		if (modeTable[i] == 1)
			continue;

		// initialize point list...
		pointCount = 0;

		// Assign window center (window centers are
		// initialized by createLattice to be the point
		// data[i])
		yk[0] = i%width;
		yk[1] = i/width;
		for(j = 0; j < N; j++)
			yk[j+2] = data[N*i+j];
		
		// Calculate the mean shift vector using the lattice
		OptLatticeMSVector(Mh, yk);
		
		// Calculate its magnitude squared
		mvAbs = 0;
		for(j = 0; j < lN; j++)
			mvAbs += Mh[j]*Mh[j];
		
		// Keep shifting window center until the magnitude squared of the
		// mean shift vector calculated at the window center location is
		// under a specified threshold (Epsilon)
		
		// NOTE: iteration count is for speed up purposes only - it
		//       does not have any theoretical importance
		iterationCount = 1;
		while((mvAbs >= EPSILON)&&(iterationCount < LIMIT))
		{
			
			// Shift window location
			for(j = 0; j < lN; j++)
				yk[j] += Mh[j];
			
			// check to see if the current mode location is in the
			// basin of attraction...

			// calculate the location of yk on the lattice
			modeCandidateX	= (int) (yk[0]+0.5);
			modeCandidateY	= (int) (yk[1]+0.5);
			modeCandidate_i	= modeCandidateY*width + modeCandidateX;

			// if mvAbs != 0 (yk did indeed move) then check
			// location basin_i in the mode table to see if
			// this data point either:
			
			// (1) has not been associated with a mode yet
			//     (modeTable[basin_i] = 0), so associate
			//     it with this one
			//
			// (2) it has been associated with a mode other
			//     than the one that this data point is converging
			//     to (modeTable[basin_i] = 1), so assign to
			//     this data point the same mode as that of basin_i

			if ((modeTable[modeCandidate_i] != 2) && (modeCandidate_i != i))
			{
				// obtain the data point at basin_i to
				// see if it is within h*TC_DIST_FACTOR of
				// of yk
				for (j = 0; j < N; j++)
					modeCandidatePoint[j] = data[N*modeCandidate_i + j];
				
				// check basin on non-spatial data spaces only
				k = 1;
				s = 0;
				diff = 0;
				while ((diff < TC_DIST_FACTOR) && (k<kp))
				{
					diff = 0;
					for (p=0; p<P[k]; p++)
					{
						el = (modeCandidatePoint[p+s]-yk[p+s+2])/h[k];
						diff += el*el;
					}
					s+=P[k];
					k++;
				}

				// if the data point at basin_i is within
				// a distance of h*TC_DIST_FACTOR of yk
				// then depending on modeTable[basin_i] perform
				// either (1) or (2)
				if (diff < TC_DIST_FACTOR)
				{
					// if the data point at basin_i has not
					// been associated to a mode then associate
					// it with the mode that this one will converge
					// to
					if (modeTable[modeCandidate_i] == 0)
					{
						// no mode associated yet so associate
						// it with this one...
						pointList[pointCount++]		= modeCandidate_i;
						modeTable[modeCandidate_i]	= 2;

					} else
					{

						// the mode has already been associated with
						// another mode, thererfore associate this one
						// mode and the modes in the point list with
						// the mode associated with data[basin_i]...

						// store the mode infor int yk using msRawData...
						for (j = 0; j < N; j++)
							yk[j+2] = msRawData[modeCandidate_i*N+j];

						// update mode table for this data point
						// indicating that a mode has been associated
						// with it
						modeTable[i] = 1;

                  // indicate that a mode has been associated
						// to this data point (data[i])
						mvAbs = -1;

						// stop mean shift calculation...
						break;
					}
				}
			}
			
			// Calculate the mean shift vector at the new
			// window location using lattice
			OptLatticeMSVector(Mh, yk);
			
			// Calculate its magnitude squared
			mvAbs = 0;
			for(j = 0; j < lN; j++)
				mvAbs += Mh[j]*Mh[j];

			// Increment interation count
			iterationCount++;
			
		}

		// if a mode was not associated with this data point
		// yet then perform a shift the window center yk one
		// last time using the mean shift vector...
		if (mvAbs >= 0)
		{
			// Shift window location
			for(j = 0; j < lN; j++)
				yk[j] += Mh[j];

         // update mode table for this data point
		   // indicating that a mode has been associated
   		// with it
			modeTable[i] = 1;
		}
		
		// associate the data point indexed by
		// the point list with the mode stored
		// by yk
		for (j = 0; j < pointCount; j++)
		{
			// obtain the point location from the
			// point list
			modeCandidate_i = pointList[j];

			// update the mode table for this point
			modeTable[modeCandidate_i] = 1;

			//store result into msRawData...
			for(k = 0; k < N; k++)
				msRawData[N*modeCandidate_i+k] = (float)(yk[k+2]);
		}


		//store result into msRawData...
		for(j = 0; j < N; j++)
			msRawData[N*i+j] = (float)(yk[j+2]);

		// Prompt user on progress
#ifdef SHOW_PROGRESS
		percent_complete = (float)(i/(float)(L))*100;
		msSys.Prompt("\r%2d%%", (int)(percent_complete + 0.5));
#endif
	
		// Check to see if the algorithm has been halted
		if((i%PROGRESS_RATE == 0)&&((ErrorStatus = msSys.Progress((float)(i/(float)(L))*(float)(0.8)))) == EL_HALT)
			break;
		
	}
	
	// Prompt user that filtering is completed
#ifdef PROMPT
#ifdef SHOW_PROGRESS
	msSys.Prompt("\r");
#endif
	msSys.Prompt("done.");
#endif
	
	// de-allocate memory
	delete [] modeCandidatePoint;
	delete [] yk;
	delete [] Mh;
	
	// done.
	return;

}

	/*/\/\/\/\/\/\/\/\/\/\/\*/
	/* Image Classification */
	/*\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Connect                                              */
/*******************************************************/
/*Classifies the regions of the mean shift filtered    */
/*image.                                               */
/*******************************************************/
/*Post:                                                */
/*      - the regions of the mean shift image have been*/
/*        classified using the private classification  */
/*        structure of the msImageProcessor Class.     */
/*        Namely, each region uniquely identified by   */
/*        its LUV color  (stored by LUV_data) and loc- */
/*        ation has been labeled and its area computed */
/*        via an eight-connected fill.                 */
/*******************************************************/

void msImageProcessor::Connect( void )
{

	//define eight connected neighbors
	neigh[0]	= 1;
	neigh[1]	= 1-width;
	neigh[2]	= -width;
	neigh[3]	= -(1+width);
	neigh[4]	= -1;
	neigh[5]	= width-1;
	neigh[6]	= width;
	neigh[7]	= width+1;

	//initialize labels and modePointCounts
	int i;
	for(i = 0; i < width*height; i++)
	{
		labels[i]			= -1;
		modePointCounts[i]	=  0;
	}

	//Traverse the image labeling each new region encountered
	int k, label = -1;
	for(i = 0; i < height*width; i++)
	{
		//if this region has not yet been labeled - label it
		if(labels[i] < 0)
		{
			//assign new label to this region
			labels[i] = ++label;    // 对每个数据点进行标识

			//copy region color into modes
			for(k = 0; k < N; k++)
            modes[(N*label)+k] = LUV_data[(N*i)+k]; // 将数据赋予给模态点，每个数据点有N维度
//			modes[(N*label)+k]	= (float)(LUV_data[(N*i)+k]);

			//populate labels with label for this specified region
			//calculating modePointCounts[label]...
			Fill(i, label); // 按照每一个点进行区域生长，类似8邻域种子填充，区域生长的方法来划分多边形
		}
	}

	//calculate region count using label
	regionCount	= label+1;

	//done.
	return;
}

/*******************************************************/
/*Fill                                                 */
/*******************************************************/
/*Given a region seed and a region label, Fill uses    */
/*the region seed to perform an eight-connected fill   */
/*for the specified region, labeling all pixels con-   */
/*tained by the region with the specified label:       */
/*label.                                               */
/*******************************************************/
/*Pre:                                                 */
/*      - regionLoc is a region seed - a pixel that is */
/*        identified as being part of the region       */
/*        labled using the label, label.               */
/*Post:                                                */
/*      - all pixels belonging to the region specified */
/*        by regionLoc (having the same integer LUV    */
/*        value specified by LUV_data) are classified  */
/*        as one region by labeling each pixel in the  */
/*        image clasification structure using label    */
/*        via an eight-connected fill.                 */
/*******************************************************/

// 通过一个种子点regionLoc，将8邻域联通内多边形内的像素点标识为同一个label，

void msImageProcessor::Fill(int regionLoc, int label)
{

	//declare variables
	int	i, k, neighLoc, neighborsFound, imageSize	= width*height;

	//Fill region starting at region location
	//using labels...

	//initialzie indexTable
	int	index		= 0;
	indexTable[0]	= regionLoc; // 存储第一个种子点

	//increment mode point counts for this region to
	//indicate that one pixel belongs to this region
	modePointCounts[label]++;

	// 如果种子填充，其8邻域的LUV数据即颜色数据与种子的颜色数据相差小于阈值，则归于种子所在的多边形
	while(true)
	{

		//assume no neighbors will be found
		neighborsFound	= 0;

		//check the eight connected neighbors at regionLoc -
		//if a pixel has similar color to that located at 
		//regionLoc then declare it as part of this region
		for(i = 0; i < 8; i++)
		{
         // no need
         /*
			//if at boundary do not check certain neighbors because
			//they do not exist...
			if((regionLoc%width == 0)&&((i == 3)||(i == 4)||(i == 5)))
				continue;
			if((regionLoc%(width-1) == 0)&&((i == 0)||(i == 1)||(i == 7)))
				continue;
         */   

			//check bounds and if neighbor has been already labeled
			neighLoc			= regionLoc + neigh[i];
			if((neighLoc >= 0)&&(neighLoc < imageSize)&&(labels[neighLoc] < 0))
			{
				for(k = 0; k < N; k++)
				{
//					if(LUV_data[(regionLoc*N)+k] != LUV_data[(neighLoc*N)+k])
               if (fabs(LUV_data[(regionLoc*N)+k]-LUV_data[(neighLoc*N)+k])>=LUV_treshold)
						break;
				}
				
				//neighbor i belongs to this region so label it and
				//place it onto the index table buffer for further
				//processing
				if(k == N)
				{
					//assign label to neighbor i
					labels[neighLoc]	= label;
					
					//increment region point count
					modePointCounts[label]++;
					
					//place index of neighbor i onto the index tabel buffer
					indexTable[++index]	= neighLoc;
					
					//indicate that a neighboring region pixel was
					//identified
					neighborsFound	= 1;
				}
			}
		}

		// 对数据进行填充
		//check the indexTable to see if there are any more
		//entries to be explored - if so explore them, otherwise
		//exit the loop - we are finished
		if(neighborsFound) // 8邻域的所有点都属于种子点所在的多边形
			regionLoc	= indexTable[index]; // 将当前点设置为下一个种子点
		else if (index > 1)
			regionLoc	= indexTable[--index];// 将前一个点设置为种子点
		else // 8邻域都不属于种子点所在的多边形
			break; //fill complete
	}
	//done.
	return;

}

	/*/\/\/\/\/\/\/\/\*/
	/*  Image Pruning */
	/*\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Build Region Adjacency Matrix                        */
/*******************************************************/
/*Constructs a region adjacency matrix.                */
/*******************************************************/
/*Pre:                                                 */
/*      - the classification data structure has been   */
/*        constructed.                                 */
/*Post:                                                */
/*      - a region adjacency matrix has been built     */
/*        using the classification data structure.     */
/*******************************************************/

// 根据labels建立RAM，也就是区域邻接表
void msImageProcessor::BuildRAM( void )
{

	//Allocate memory for region adjacency matrix if it hasn't already been allocated
	if((!raList)&&((!(raList = new RAList [regionCount]))||(!(raPool = new RAList [NODE_MULTIPLE*regionCount]))))
	{
		ErrorHandler("msImageProcessor", "Allocate", "Not enough memory.");
		return;
	}

	//initialize the region adjacency list
	int i;
	// ralist是一个邻接表，每个结点有如下的属性:edgeStrength,edgePixelCount,label,next
	for(i = 0; i < regionCount; i++)
	{
		raList[i].edgeStrength		= 0;
		raList[i].edgePixelCount	= 0;
		raList[i].label				= i; // 这个label与那个region的label有何区别
		raList[i].next				= NULL;
	}

	//initialize RAM free list
	freeRAList	= raPool; // 这个不是邻接表，是单链表
	for(i = 0; i < NODE_MULTIPLE*regionCount-1; i++)
	{
		raPool[i].edgeStrength		= 0;  // 边缘强度
		raPool[i].edgePixelCount	= 0;  // 边缘像素个数
		raPool[i].next = &raPool[i+1];    // 
	}
	raPool[NODE_MULTIPLE*regionCount-1].next	= NULL;

	//traverse the labeled image building
	//the RAM by looking to the right of
	//and below the current pixel location thus
	//determining if a given region is adjacent
	//to another
	int	 curLabel, rightLabel, bottomLabel, exists;
	RAList	*raNode1, *raNode2, *oldRAFreeList;


	/************************************************************************/
	/* 注释掉的代码--开头                                                   */
	/************************************************************************/
	/*
	int j;
	for(i = 0; i < height - 1; i++)
	{
		//check the right and below neighbors
		//for pixel locations whose x < width - 1
		for(j = 0; j < width - 1; j++)
		{
			//calculate pixel labels
			curLabel	= labels[i*width+j    ];	//current pixel
			rightLabel	= labels[i*width+j+1  ];	//right   pixel
			bottomLabel	= labels[(i+1)*width+j];	//bottom  pixel

			//check to the right, if the label of
			//the right pixel is not the same as that
			//of the current one then region[j] and region[j+1]
			//are adjacent to one another - update the RAM
			if(curLabel != rightLabel)
			{
				//obtain RAList object from region adjacency free
				//list
				raNode1			= freeRAList; // 
				raNode2			= freeRAList->next;

				//keep a pointer to the old region adj. free
				//list just in case nodes already exist in respective
				//region lists
				oldRAFreeList	= freeRAList;

				//update region adjacency free list
				freeRAList		= freeRAList->next->next;

				//populate RAList nodes
				raNode1->label	= curLabel; // 当前多边形
				raNode2->label	= rightLabel;// 右边多边形

				//insert nodes into the RAM
				exists			= 0;
				raList[curLabel  ].Insert(raNode2);// 右边多边形与当前多边形相邻，插入RAM，也就是插入图的邻接表
				exists			= raList[rightLabel].Insert(raNode1);

				//if the node already exists then place
				//nodes back onto the region adjacency
				//free list
				if(exists)
					freeRAList = oldRAFreeList; 

			}

			//check below, if the label of
			//the bottom pixel is not the same as that
			//of the current one then region[j] and region[j+width]
			//are adjacent to one another - update the RAM
			if(curLabel != bottomLabel)
			{
				//obtain RAList object from region adjacency free
				//list
				raNode1			= freeRAList;
				raNode2			= freeRAList->next;

				//keep a pointer to the old region adj. free
				//list just in case nodes already exist in respective
				//region lists
				oldRAFreeList	= freeRAList;

				//update region adjacency free list
				freeRAList		= freeRAList->next->next;

				//populate RAList nodes
				raNode1->label	= curLabel; // 计算当前标识
				raNode2->label	= bottomLabel; // 计算底部标识

				//insert nodes into the RAM
				exists			= 0;
				raList[curLabel].Insert(raNode2);
				exists			= raList[bottomLabel].Insert(raNode1);

				//if the node already exists then place
				//nodes back onto the region adjacency
				//free list
				if(exists)
					freeRAList = oldRAFreeList;

			}

		}// for j为0到width-2



		//check only to the bottom neighbors of the right boundary
		//pixels...

		//calculate pixel locations (j = width-1)，最右边一列，只能往下
		curLabel	= labels[i*width+j    ];	//current pixel
		bottomLabel = labels[(i+1)*width+j];	//bottom  pixel

		//check below, if the label of
		//the bottom pixel is not the same as that
		//of the current one then region[j] and region[j+width]
		//are adjacent to one another - update the RAM
		if(curLabel != bottomLabel)
		{
			//obtain RAList object from region adjacency free
			//list
			raNode1			= freeRAList;
			raNode2			= freeRAList->next;
			
			//keep a pointer to the old region adj. free
			//list just in case nodes already exist in respective
			//region lists
			oldRAFreeList	= freeRAList;
			
			//update region adjacency free list
			freeRAList		= freeRAList->next->next;
			
			//populate RAList nodes
			raNode1->label	= curLabel;
			raNode2->label	= bottomLabel;
			
			//insert nodes into the RAM
			exists			= 0;
			raList[curLabel  ].Insert(raNode2);
			exists			= raList[bottomLabel].Insert(raNode1);
			
			//if the node already exists then place
			//nodes back onto the region adjacency
			//free list
			if(exists)
				freeRAList = oldRAFreeList;

		}
	}// for i = 0到height-2

	//check only to the right neighbors of the bottom boundary
	//pixels...

	//check the right for pixel locations whose x < width - 1// 最下面一行，只能往右
	for(j = 0; j < width - 1; j++)
	{
		//calculate pixel labels (i = height-1)
		curLabel	= labels[i*width+j    ];	//current pixel
		rightLabel	= labels[i*width+j+1  ];	//right   pixel
		
		//check to the right, if the label of
		//the right pixel is not the same as that
		//of the current one then region[j] and region[j+1]
		//are adjacent to one another - update the RAM
		if(curLabel != rightLabel)
		{
			//obtain RAList object from region adjacency free
			//list
			raNode1			= freeRAList;
			raNode2			= freeRAList->next;

			//keep a pointer to the old region adj. free
			//list just in case nodes already exist in respective
			//region lists
			oldRAFreeList	= freeRAList;
			
			//update region adjacency free list
			freeRAList		= freeRAList->next->next;
			
			//populate RAList nodes
			raNode1->label	= curLabel;
			raNode2->label	= rightLabel;
			
			//insert nodes into the RAM
			exists			= 0;
			raList[curLabel  ].Insert(raNode2);
			exists			= raList[rightLabel].Insert(raNode1);
			
			//if the node already exists then place
			//nodes back onto the region adjacency
			//free list
			if(exists)
				freeRAList = oldRAFreeList;

		}

	}

/************************************************************************/
/* 注释掉的代码--结束                                                   */
/************************************************************************/

if (!regionList)
{
	regionList = GetBoundariesAndRegions();
}

int EdgePtNums = regionList->GetTotalBoundaryPtCount();

int pos, xPos, yPos; // 获取多边形的标识
int *boundaryPtIdx = regionList->GetBoundaryIndeces(0);
curLabel = -1;
rightLabel = -1; 
bottomLabel = -1;
for (i = 0; i < EdgePtNums; i ++)
{
	pos = boundaryPtIdx[i];
	xPos = pos % width; // 边界点的列号
	yPos = pos / width; // 边界点的行号
	curLabel = labels[pos];
	if (xPos != width - 1)
	{
		rightLabel = labels[pos+1];
	}
	if (yPos != height - 1)
	{
		bottomLabel = labels[pos+width];
	}
	if(curLabel != rightLabel && rightLabel != -1)
	{
		//obtain RAList object from region adjacency free
		//list
		raNode1			= freeRAList; // 
		raNode2			= freeRAList->next;
		
		//keep a pointer to the old region adj. free
		//list just in case nodes already exist in respective
		//region lists
		oldRAFreeList	= freeRAList;
		
		//update region adjacency free list
		freeRAList		= freeRAList->next->next;
		
		//populate RAList nodes
		raNode1->label	= curLabel; // 当前多边形
		raNode2->label	= rightLabel;// 右边多边形
		
		//insert nodes into the RAM
		exists			= 0;
		raList[curLabel  ].Insert(raNode2);// 右边多边形与当前多边形相邻，插入RAM，也就是插入图的邻接表
		exists			= raList[rightLabel].Insert(raNode1);
		
		//if the node already exists then place
		//nodes back onto the region adjacency
		//free list
		if(exists)
			freeRAList = oldRAFreeList; 
		
	}
	
	//check below, if the label of
	//the bottom pixel is not the same as that
	//of the current one then region[j] and region[j+width]
	//are adjacent to one another - update the RAM
	if(curLabel != bottomLabel &&  bottomLabel != -1)
	{
		//obtain RAList object from region adjacency free
		//list
		raNode1			= freeRAList;
		raNode2			= freeRAList->next;
		
		//keep a pointer to the old region adj. free
		//list just in case nodes already exist in respective
		//region lists
		oldRAFreeList	= freeRAList;
		
		//update region adjacency free list
		freeRAList		= freeRAList->next->next;
		
		//populate RAList nodes
		raNode1->label	= curLabel; // 计算当前标识
		raNode2->label	= bottomLabel; // 计算底部标识
		
		//insert nodes into the RAM
		exists			= 0;
		raList[curLabel].Insert(raNode2);
		exists			= raList[bottomLabel].Insert(raNode1);
		
		//if the node already exists then place
		//nodes back onto the region adjacency
		//free list
		if(exists)
			freeRAList = oldRAFreeList;
		
		}

}


	//done.
	return;

}

/*******************************************************/
/*Destroy Region Adjacency Matrix                      */
/*******************************************************/
/*Destroy a region adjacency matrix.                   */
/*******************************************************/
/*Post:                                                */
/*      - the region adjacency matrix has been destr-  */
/*        oyed: (1) its memory has been de-allocated,  */
/*        (2) the RAM structure has been initialize    */
/*        for re-use.                                  */
/*******************************************************/

void msImageProcessor::DestroyRAM( void )
{

	//de-allocate memory for region adjaceny list
	if (raList)				delete [] raList; 
	if (raPool)				delete [] raPool; 

	//initialize region adjacency matrix
	raList				= NULL;
	freeRAList			= NULL;
	raPool				= NULL;

	//done.
	return;

}

/*******************************************************/
/*Transitive Closure                                   */
/*******************************************************/
/*Applies transitive closure to the RAM updating       */
/*labels, modes and modePointCounts to reflect the new */
/*set of merged regions resulting from transitive clo- */
/*sure.                                                */
/*******************************************************/
/*Post:                                                */
/*      - transitive closure has been applied to the   */
/*        regions classified by the RAM and labels,    */
/*        modes and modePointCounts have been updated  */
/*        to reflect the new set of mergd regions res- */
/*        ulting from transitive closure.              */
/*******************************************************/
// 什么是RAM，其是如何应用起来的，region adjacent matrix

void msImageProcessor::TransitiveClosure( void )
{

	//Step (1):

	// Build RAM using classifiction structure originally
	// generated by the method GridTable::Connect()
	BuildRAM();

	//Step (1a):
	//Compute weights of weight graph using confidence map
	//(if defined)
	if(weightMapDefined)	ComputeEdgeStrengths(); // 计算边缘强度

	//Step (2):

	//Treat each region Ri as a disjoint set:

	// - attempt to join Ri and Rj for all i != j that are neighbors and
	//   whose associated modes are a normalized distance of < 0.5 from one
	//   another

	// - the label of each region in the raList is treated as a pointer to the
	//   canonical element of that region (e.g. raList[i], initially has raList[i].label = i,
	//   namely each region is initialized to have itself as its canonical element).



	// 如果两个多边形相邻且她们InWindow（即whose associated modes are a normalized distance of < 0.5），则将其合并

	// ************************************************************************遍历RAM，尝试将raList[i]与其相邻的多边形进行合并
	//Traverse RAM attempting to join raList[i] with its neighbors...
	int		i, iCanEl, neighCanEl;
	float	threshold;   // 阈值
	RAList	*neighbor;   // 邻接多边形
	for(i = 0; i < regionCount; i++)
	{
		//acquire first neighbor in region adjacency list pointed to
		//by raList[i]
		neighbor	= raList[i].next;

		//compute edge strength threshold using global and local
		//epsilon
		if(epsilon > raList[i].edgeStrength)
			threshold   = epsilon;
		else
			threshold   = raList[i].edgeStrength;

		//traverse region adjacency list of region i, attempting to join
		//it with regions whose mode is a normalized distance < 0.5 from
		//that of region i...
		// 依次遍历其邻居，寻找合适的进行合并
		while(neighbor)
		{
			//attempt to join region and neighbor...
			if((InWindow(i, neighbor->label))&&(neighbor->edgeStrength < epsilon))
			{
				//region i and neighbor belong together so join them
				//by:

				// (1) find the canonical element of region i
				iCanEl		= i;
				while(raList[iCanEl].label != iCanEl)
					iCanEl		= raList[iCanEl].label;

				// (2) find the canonical element of neighboring region
				neighCanEl	= neighbor->label;
				while(raList[neighCanEl].label != neighCanEl)
					neighCanEl	= raList[neighCanEl].label;  // zc的那个rank传递其实就是一个传递闭包来进行合并相邻多边形的过程

				// if the canonical elements of are not the same then assign
				// the canonical element having the smaller label to be the parent
				// of the other region...
				if(iCanEl < neighCanEl)
					raList[neighCanEl].label	= iCanEl; // 将neighCanEl归并到iCanEl？将大的编号归并到小的编号中？
				else
				{
					//must replace the canonical element of previous
					//parent as well
					raList[raList[iCanEl].label].label	= neighCanEl;// 将iCanEl之流归并到neighCanEl？

					//re-assign canonical element
					raList[iCanEl].label				= neighCanEl;
				}
			}

			//check the next neighbor...
			neighbor	= neighbor->next; //找下一个邻居

		}
	}

	// Step (3):

	// Level binary trees formed by canonical elements
	// 通过canonical elements形成层次二叉树
	for(i = 0; i < regionCount; i++)
	{
		iCanEl	= i;
		while(raList[iCanEl].label != iCanEl)
			iCanEl	= raList[iCanEl].label; // 自下往上，也就是自叶子往根部回溯
		raList[i].label	= iCanEl;
	}

	// Step (4):

	//Traverse joint sets, relabeling image. // 遍历各区域，重新label

	// (a)

	// Accumulate modes and re-compute point counts using canonical
	// elements generated by step 2.

	//allocate memory for mode and point count temporary buffers...
	float	*modes_buffer	= new float	[N*regionCount];
	int		*MPC_buffer		= new int	[regionCount];

	//initialize buffers to zero
	for(i = 0; i < regionCount; i++)
		MPC_buffer[i]	= 0;
	for(i = 0; i < N*regionCount; i++)
		modes_buffer[i]	= 0;

	//traverse raList accumulating modes and point counts
	//using canoncial element information...
	int k, iMPC;
	for(i = 0; i < regionCount; i++)
	{

		//obtain canonical element of region i
		iCanEl	= raList[i].label;

		//obtain mode point count of region i
		iMPC	= modePointCounts[i];

		//accumulate modes_buffer[iCanEl]
		for(k = 0; k < N; k++)
			modes_buffer[(N*iCanEl)+k] += iMPC*modes[(N*i)+k]; // modes类似于data的均值，i多边形N维度的均值

		//accumulate MPC_buffer[iCanEl]
		MPC_buffer[iCanEl] += iMPC;

	}

	// (b)

	// Re-label new regions of the image using the canonical
	// element information generated by step (2)

	// Also use this information to compute the modes of the newly
	// defined regions, and to assign new region point counts in
	// a consecute manner to the modePointCounts array

	//allocate memory for label buffer
	int	*label_buffer	= new int [regionCount];

	//initialize label buffer to -1
	for(i = 0; i < regionCount; i++)
		label_buffer[i]	= -1;

	//traverse raList re-labeling the regions
	int	label = -1;
	for(i = 0; i < regionCount; i++)
	{
		//obtain canonical element of region i
		iCanEl	= raList[i].label;
		if(label_buffer[iCanEl] < 0)
		{
			//assign a label to the new region indicated by canonical
			//element of i
			label_buffer[iCanEl]	= ++label; // 同一个iCanEl，其实就是同一个进行合并的

			//recompute mode storing the result in modes[label]...
			iMPC	= MPC_buffer[iCanEl];
			for(k = 0; k < N; k++)
				modes[(N*label)+k]	= (modes_buffer[(N*iCanEl)+k])/(iMPC); // 均值，可恶的均值，更新多边形合并后的各多边形的均值modes

			//assign a corresponding mode point count for this region into
			//the mode point counts array using the MPC buffer...
			modePointCounts[label]	= MPC_buffer[iCanEl];
		}
	}

	//re-assign region count using label counter
	int	oldRegionCount	= regionCount;
	regionCount	= label+1;

	// (c)

	// Use the label buffer to reconstruct the label map, which specified
	// the new image given its new regions calculated above

	/************************************************************labels的华丽转变，众里寻你千百度，o(∩_∩)o 哈哈**********************************************************************/	
	for(i = 0; i < height*width; i++)
		labels[i]	= label_buffer[raList[labels[i]].label]; // 合并后更新label

	//de-allocate memory
	delete [] modes_buffer;
	delete [] MPC_buffer;
	delete [] label_buffer;

	//done.
	return;

}

/*******************************************************/
/*Compute Edge Strengths                               */
/*******************************************************/
/*Computes the a weight for each link in the region    */
/*graph maintined by the RAM, resulting in a weighted  */
/*graph in which the weights consist of a confidence   */
/*between zero and one indicating if the regions are   */
/*separated by a strong or weak edge.                  */
/*******************************************************/
/*Post:                                                */
/*      - an edge strength has been computed between   */
/*        each region of the image and placed as a     */
/*        weight in the RAM to be used during transi-  */
/*        tive closure.                                */
/*******************************************************/

void msImageProcessor::ComputeEdgeStrengths( void )
{

	//initialize visit table - used to keep track
	//of which pixels have already been visited such
	//as not to contribute their strength value to
	//a boundary sum multiple times...
	memset(visitTable, 0, L*sizeof(unsigned char));

	//traverse labeled image computing edge strengths
	//(excluding image boundary)...
	int    x, dp, curLabel, rightLabel, bottomLabel;
	RAList *curRegion;
	/************************************************************************/
	/* 注释掉一下代码，对代码进行优化，以提高速度 ,这代码与下面的代码意义一样                          */
	/************************************************************************/
/*
int y;
	for(y = 1; y < height-1; y++)
	{
		for(x = 1; x < width-1; x++)
		{
			//compute data point location using x and y
			dp = y*width + x;

			//obtain labels at different pixel locations
			curLabel	= labels[dp      ];	//current pixel
			rightLabel	= labels[dp+1    ];	//right   pixel
			bottomLabel	= labels[dp+width];	//bottom  pixel

			//check right and bottom neighbor to see if there is a
			//change in label then we are at an edge therefore record
			//the edge strength at this edge accumulating its value
			//in the RAM...
			// 如果扫描RegionList中的indexTable多好了，存储了所有多边形的边界像素的索引，我们只要比较边界像素的右边以及往下就ok了
			if(curLabel != rightLabel) // 这些其实都是边缘像素，在buildRAM也是的，邻接表不是建立好了吗？为什么不从扫描邻接表开始啊，邻接表的label没有变啊！！因为需要weightMap
			{
				// 遍历RAM，获得右边的Region
				//traverse into RAM...
				curRegion = &raList[curLabel]; // 先计算邻接表每一个链条后面结点的edgeStrength属性
				while((curRegion)&&(curRegion->label != rightLabel))
					curRegion = curRegion->next;

				//this should not occur...
				assert(curRegion);

				//accumulate edge strength，为什么此处只计算rightlabel或是bottomlabel的Region的edgeStrength，难道是计算单次？
				curRegion->edgeStrength   += weightMap[dp] + weightMap[dp+1]; // 右相邻多边形的edgeStrength计算
				curRegion->edgePixelCount += 2;
			}

			if(curLabel != bottomLabel)// 这些其实都是边缘像素，在buildRAM也是的，难道是最后才确定region，之前都用对所有点label标识
			{
				//traverse into RAM...
				curRegion = &raList[curLabel];
				// 遍历RAM，获得下面相邻的Region
				while((curRegion)&&(curRegion->label != bottomLabel))
					curRegion = curRegion->next;

				//this should not occur...
				assert(curRegion);

				//accumulate edge strength
				if(curLabel == rightLabel)
				{
					curRegion->edgeStrength   += weightMap[dp] + weightMap[dp+width];
					curRegion->edgePixelCount += 2;
				} 
				else
				{
					curRegion->edgeStrength	  += weightMap[dp+width];
					curRegion->edgePixelCount += 1;
				}

			}
		}
	}
	*/
	/************************************************************************/
	/* 注释掉一下代码，对代码进行优化，以提高速度                           */
	/************************************************************************/


	// 以边界的数据进行处理，问题是要获取边界，并将边界进行存储，每一次进行计算RAM都要进行一次获取边界，因为lables发生改变
	// 我们假设在buildRAM的时候已经进行了GetBoundaries。也就是获得regionList通过边界对其进行处理

	/************************************************************************/
	/* *************************自定义代码-start*****************************/
	/************************************************************************/
	if (!regionList)
	{
		regionList = GetBoundariesAndRegions();
	}

	int EdgePtNums = regionList->GetTotalBoundaryPtCount();

	int xPos, yPos; // 获取多边形的标识
	int *boundaryPtIdx = regionList->GetBoundaryIndeces(0);

	for (int i = 0; i < EdgePtNums; i ++)
	{
		dp = boundaryPtIdx[i]; // 获取边界点
		
		xPos = dp % width; // 边界点的列号
		yPos = dp / width; // 边界点的行号
		if (xPos > 0 && xPos < width - 1 && yPos > 0 && yPos < height - 1)
		{
			curLabel = labels[dp];
			rightLabel = labels[dp+1];
			bottomLabel = labels[dp+width];
			//check right and bottom neighbor to see if there is a
			//change in label then we are at an edge therefore record
			//the edge strength at this edge accumulating its value
			//in the RAM...
			// 如果扫描RegionList中的indexTable多好了，存储了所有多边形的边界像素的索引，我们只要比较边界像素的右边以及往下就ok了
			if(curLabel != rightLabel) // 这些其实都是边缘像素，在buildRAM也是的，邻接表不是建立好了吗？为什么不从扫描邻接表开始啊，邻接表的label没有变啊！！因为需要weightMap
			{
				// 遍历RAM，获得右边的Region
				//traverse into RAM...
				curRegion = &raList[curLabel]; // 先计算邻接表每一个链条后面结点的edgeStrength属性
				while((curRegion)&&(curRegion->label != rightLabel))
					curRegion = curRegion->next;
				
				//this should not occur...
				assert(curRegion);
				
				//accumulate edge strength，为什么此处只计算rightlabel或是bottomlabel的Region的edgeStrength，难道是计算单次？
				curRegion->edgeStrength   += weightMap[dp] + weightMap[dp+1]; // 右相邻多边形的edgeStrength计算
				curRegion->edgePixelCount += 2;
			}
			
			if(curLabel != bottomLabel)// 这些其实都是边缘像素，在buildRAM也是的，难道是最后才确定region，之前都用对所有点label标识
			{
				//traverse into RAM...
				curRegion = &raList[curLabel];
				// 遍历RAM，获得下面相邻的Region
				while((curRegion)&&(curRegion->label != bottomLabel))
					curRegion = curRegion->next;
				
				//this should not occur...
				assert(curRegion);
				
				//accumulate edge strength
				if(curLabel == rightLabel)
				{
					curRegion->edgeStrength   += weightMap[dp] + weightMap[dp+width];
					curRegion->edgePixelCount += 2;
				} 
				else
				{
					curRegion->edgeStrength	  += weightMap[dp+width];
					curRegion->edgePixelCount += 1;
				}
				
			}
		}
	}
	/************************************************************************/
	/* *************************自定义代码-end******************************/
	/************************************************************************/


	//compute strengths using accumulated strengths obtained above...
	RAList *neighborRegion;
	float	edgeStrength;
	int		edgePixelCount;
	// 先计算与raList[x]相邻的多边形与其相邻的那一部分的edgeStrength
	// 再统计raList[x]所有相邻多边形的edgeStrength，求平均值来得到其edgeStrength
	// 先部分后整体，先周边再统筹
	for(x = 0; x < regionCount; x++)
	{
		//traverse the region list of the current region
		curRegion	= &raList[x];
		curRegion	= curRegion->next;
		while(curRegion)
		{
			//with the assumption that regions having a smaller
			//label in the current region list have already
			//had their edge strengths computed, only compute
			//edge strengths for the regions whose label is greater
			//than x, the current region (region list) under
			//consideration...
			// 我们假设其label比x小的邻接多边形的edgeStrength都计算好了，我们从其label比x大的邻接多边形计算
			curLabel = curRegion->label;
			if(curLabel > x)
			{
				//obtain pointer to the element identifying the
				//current region in the neighbors region list...
				neighborRegion = &raList[curLabel];// 从其比x大的邻接多边形来了解x的edgeStrength
				while((neighborRegion)&&(neighborRegion->label != x))
					neighborRegion = neighborRegion->next; // 
				
				//this should not occur...
				assert(neighborRegion);
				
				//compute edge strengths using accumulated confidence
				//value and pixel count
				if((edgePixelCount = curRegion->edgePixelCount + neighborRegion->edgePixelCount) != 0)
				{
					//compute edge strength
					edgeStrength	= curRegion->edgeStrength + neighborRegion->edgeStrength;
					edgeStrength	/= edgePixelCount;
					
					//store edge strength and pixel count for corresponding regions
					curRegion->edgeStrength		= neighborRegion->edgeStrength		= edgeStrength;
					curRegion->edgePixelCount	= neighborRegion->edgePixelCount	= edgePixelCount;
				}
			}

			//traverse to the next region in the region adjacency list
			//of the current region x
			curRegion = curRegion->next;

		}
	}

	//compute average edge strength amongst the edges connecting
	//it to each of its neighbors
	// 计算平均的边缘强度
	int numNeighbors;
	for(x = 0; x < regionCount; x++)
	{
		//traverse the region list of the current region
		//accumulating weights
		curRegion		= &raList[x];      // 
		curRegion		= curRegion->next; //
		edgeStrength	= 0;
		numNeighbors	= 0;

		while(curRegion)
		{
			numNeighbors++; // 相邻多边形的个数
			edgeStrength   += curRegion->edgeStrength; // 边缘强度
			curRegion		= curRegion->next; // 当前的多边形
		}

		//divide by the number of regions connected
		//to the current region
		if(numNeighbors) edgeStrength /= numNeighbors;  // 计算平均边缘强度=相邻多边形的边缘强度/相邻多边形的个数

		//store the result in the raList for region
		//x
		raList[x].edgeStrength = edgeStrength;
	}


	//traverse raList and output the resulting list
	//to a file

	//done.
	return;

}

/*******************************************************/
/*Prune                                                */
/*******************************************************/
/*Prunes regions from the image whose pixel density    */
/*is less than a specified threshold.                  */
/*******************************************************/
/*Pre:                                                 */
/*      - minRegion is the minimum allowable pixel de- */
/*        nsity a region may have without being pruned */
/*        from the image                               */
/*Post:                                                */
/*      - regions whose pixel density is less than     */
/*        or equal to minRegion have been pruned from  */
/*        the image.                                   */
/*******************************************************/

void msImageProcessor::Prune(int minRegion)
{
	
	//Allocate Memory for temporary buffers...
	
	//allocate memory for mode and point count temporary buffers...
	float	*modes_buffer	= new float	[N*regionCount];
	int		*MPC_buffer		= new int	[regionCount];
	
	//allocate memory for label buffer
	int	*label_buffer		= new int	[regionCount];
	
	//Declare variables
	int		i, k, candidate, iCanEl, neighCanEl, iMPC, label, oldRegionCount, minRegionCount;
	double	minSqDistance, neighborDistance;
	RAList	*neighbor;
	
	//Apply pruning algorithm to classification structure, removing all regions whose area
	//is under the threshold area minRegion (pixels)
	do
	{
		//Assume that no region has area under threshold area  of 
		minRegionCount	= 0;		

		//Step (1):
		
		// Build RAM using classifiction structure originally
		// generated by the method GridTable::Connect()
		BuildRAM();
		
		// Step (2):
		
		// Traverse the RAM joining regions whose area is less than minRegion (pixels)
		// with its respective candidate region.
		
		// A candidate region is a region that displays the following properties:
		
		//	- it is adjacent to the region being pruned
		
		//  - the distance of its mode is a minimum to that of the region being pruned
		//    such that or it is the only adjacent region having an area greater than
		//    minRegion
		
		for(i = 0; i < regionCount; i++)
		{
			//if the area of the ith region is less than minRegion
			//join it with its candidate region...

			//*******************************************************************************

			//Note: Adjust this if statement if a more sophisticated pruning criterion
			//      is desired. Basically in this step a region whose area is less than
			//      minRegion is pruned by joining it with its "closest" neighbor (in color).
			//      Therefore, by placing a different criterion for fusing a region the
			//      pruning method may be altered to implement a more sophisticated algorithm.

			//*******************************************************************************

			if(modePointCounts[i] < minRegion)
			{
				//update minRegionCount to indicate that a region
				//having area less than minRegion was found
				minRegionCount++;

				//obtain a pointer to the first region in the
				//region adjacency list of the ith region...
				neighbor	= raList[i].next;
				
				//calculate the distance between the mode of the ith
				//region and that of the neighboring region...
				candidate		= neighbor->label;
				minSqDistance	= SqDistance(i, candidate);
				
				//traverse region adjacency list of region i and select
				//a candidate region
				neighbor	= neighbor->next;
				// 从邻居中找到最近的进行合并，大的label合并到小的label中
				while(neighbor)
				{

					//calculate the square distance between region i
					//and current neighbor...
					neighborDistance = SqDistance(i, neighbor->label);

					//if this neighbors square distance to region i is less
					//than minSqDistance, then select this neighbor as the
					//candidate region for region i
					if(neighborDistance < minSqDistance)
					{
						minSqDistance	= neighborDistance;
						candidate		= neighbor->label;
					}

					//traverse region list of region i
					neighbor	= neighbor->next;

				}

/****************************************************join与relable都可以自己成一个函数，这样可以Dont Repeat Yourself************************************/
				//join region i with its candidate region:

				// (1) find the canonical element of region i
				iCanEl		= i;
				while(raList[iCanEl].label != iCanEl)
					iCanEl		= raList[iCanEl].label;

				// (2) find the canonical element of neighboring region
				neighCanEl	= candidate;
				while(raList[neighCanEl].label != neighCanEl)
					neighCanEl	= raList[neighCanEl].label;

				// if the canonical elements of are not the same then assign
				// the canonical element having the smaller label to be the parent
				// of the other region...
				if(iCanEl < neighCanEl)
					raList[neighCanEl].label	= iCanEl;
				else
				{
					//must replace the canonical element of previous
					//parent as well
					raList[raList[iCanEl].label].label	= neighCanEl;

					//re-assign canonical element
					raList[iCanEl].label				= neighCanEl;
				}
			}// if
		}// for

		// Step (3):
		
		// Level binary trees formed by canonical elements
		for(i = 0; i < regionCount; i++)
		{
			iCanEl	= i;
			while(raList[iCanEl].label != iCanEl)
				iCanEl	= raList[iCanEl].label;
			raList[i].label	= iCanEl;
		}
		
		// Step (4):
		
		//Traverse joint sets, relabeling image.
		
		// Accumulate modes and re-compute point counts using canonical
		// elements generated by step 2.
		
		//initialize buffers to zero
		for(i = 0; i < regionCount; i++)
			MPC_buffer[i]	= 0;
		for(i = 0; i < N*regionCount; i++)
			modes_buffer[i]	= 0;
		
		//traverse raList accumulating modes and point counts
		//using canoncial element information...
		for(i = 0; i < regionCount; i++)
		{
			
			//obtain canonical element of region i
			iCanEl	= raList[i].label;
			
			//obtain mode point count of region i
			iMPC	= modePointCounts[i];
			
			//accumulate modes_buffer[iCanEl]
			for(k = 0; k < N; k++)
				modes_buffer[(N*iCanEl)+k] += iMPC*modes[(N*i)+k];
			
			//accumulate MPC_buffer[iCanEl]
			MPC_buffer[iCanEl] += iMPC;
			
		}
		
		// (b)
		
		// Re-label new regions of the image using the canonical
		// element information generated by step (2)
		
		// Also use this information to compute the modes of the newly
		// defined regions, and to assign new region point counts in
		// a consecute manner to the modePointCounts array
		
		//initialize label buffer to -1
		for(i = 0; i < regionCount; i++)
			label_buffer[i]	= -1;
		
		//traverse raList re-labeling the regions
		label = -1;
		for(i = 0; i < regionCount; i++)
		{
			//obtain canonical element of region i
			iCanEl	= raList[i].label;
			if(label_buffer[iCanEl] < 0)
			{
				//assign a label to the new region indicated by canonical
				//element of i
				label_buffer[iCanEl]	= ++label;
				
				//recompute mode storing the result in modes[label]...
				iMPC	= MPC_buffer[iCanEl];
				for(k = 0; k < N; k++)
					modes[(N*label)+k]	= (modes_buffer[(N*iCanEl)+k])/(iMPC);
				
				//assign a corresponding mode point count for this region into
				//the mode point counts array using the MPC buffer...
				modePointCounts[label]	= MPC_buffer[iCanEl];
			}
		}
		
		//re-assign region count using label counter
		oldRegionCount	= regionCount;
		regionCount		= label+1;
		
		// (c)
		
		// Use the label buffer to reconstruct the label map, which specified
		// the new image given its new regions calculated above
		
		for(i = 0; i < height*width; i++)
			labels[i]	= label_buffer[raList[labels[i]].label];

		
	}	while(minRegionCount > 0);  // do……while

	//de-allocate memory
	delete [] modes_buffer;
	delete [] MPC_buffer;
	delete [] label_buffer;
	
	//done.
	return;
	
}

/*******************************************************/
/*Define Boundaries                                    */
/*******************************************************/
/*Defines the boundaries for each region of the segm-  */
/*ented image storing the result into a region list    */
/*object.                                              */
/*******************************************************/
/*Pre:                                                 */
/*      - the image has been segmented and a classifi- */
/*        cation structure has been created for this   */
/*        image                                        */
/*Post:                                                */
/*      - the boundaries of the segmented image have   */
/*        been defined and the boundaries of each reg- */
/*        ion has been stored into a region list obj-  */
/*        ect.                                         */
/*******************************************************/

void msImageProcessor::DefineRegions( void )
{

	//declare and allocate memory for boundary map and count
	int	*boundaryMap,	*boundaryCount;
	if((!(boundaryMap = new int [L]))||(!(boundaryCount = new int [regionCount])))
		ErrorHandler("msImageProcessor", "DefineRegions", "Not enough memory.");

	//initialize boundary map and count
	int i;
	for(i = 0; i < L; i++)
		boundaryMap[i]		= -1;  // 用一个boundaryMap来标识boundary
	for(i = 0; i < regionCount; i++)
		boundaryCount[i]	=  0; // 记录每一个region的boundary的点数

	//initialize and declare total boundary count -
	//the total number of boundary pixels present in
	//the segmented image
	int	totalBoundaryCount	= 0; //分割图像的整个边界的像素的个数

	//traverse the image checking the right and bottom
	//four connected neighbors of each pixel marking
	//boundary map with the boundaries of each region and
	//incrementing boundaryCount using the label information

	//***********************************************************************
	//***********************************************************************

	// 扫描图像，从上往下，从左至右。每个像素的4连通区域，每个多边形的边界组成的边界map，
	// 用label信息来增加boundaryCount

	int		j, label, dataPoint;

	//first row (every pixel is a boundary pixel)
	//第一行
	for(i = 0; i < width; i++)
	{
			boundaryMap[i]		= label	= labels[i];
			boundaryCount[label]++; // label是多边形的标识。强大的labels，比indexTable更强大
			totalBoundaryCount++;
	}

	//define boundaries for all rows except for the first
	//and last one...
	for(i = 1; i < height - 1; i++)
	{
		//每行最左边的为边界，予以标识
		//mark the first pixel in an image row as an image boundary...
		dataPoint				= i*width;
		boundaryMap[dataPoint]	= label	= labels[dataPoint];
		boundaryCount[label]++;
		totalBoundaryCount++;

		for(j = 1; j < width - 1; j++)
		{
			//define datapoint and its right and bottom
			//four connected neighbors
			dataPoint		= i*width+j;

			//check four connected neighbors if they are
			//different this pixel is a boundary pixel
			label	= labels[dataPoint]; // 检查该点的四邻域如果该点的标识不等于其四领域中的一个，则其为边界点
			if((label != labels[dataPoint-1])    ||(label != labels[dataPoint+1])||
			   (label != labels[dataPoint-width])||(label != labels[dataPoint+width]))
			{
				boundaryMap[dataPoint]		= label	= labels[dataPoint];
				boundaryCount[label]++;
				totalBoundaryCount++;
			}
		}

		//mark the last pixel in an image row as an image boundary...
		dataPoint				= (i+1)*width-1;
		boundaryMap[dataPoint]	= label	= labels[dataPoint];
		boundaryCount[label]++;
		totalBoundaryCount++;

	}

	//last row (every pixel is a boundary pixel) (i = height-1)
	register int	start	= (height-1)*width, stop = height*width;
	for(i = start; i < stop; i++)
	{
		boundaryMap[i]		= label	= labels[i];
		boundaryCount[label]++;
		totalBoundaryCount++;
	}


	int	*boundaryBuffer	= new int [totalBoundaryCount], *boundaryIndex	= new int [regionCount];

	//use boundary count to initialize boundary index...
	int counter = 0;
	for(i = 0; i < regionCount; i++)
	{
		boundaryIndex[i]	= counter;
		counter			   += boundaryCount[i];
	}

	//traverse boundary map placing the boundary pixel
	//locations into the boundaryBuffer
	for(i = 0; i < L; i++)
	{
		//if its a boundary pixel store it into
		//the boundary buffer
		if((label = boundaryMap[i]) >= 0)
		{
			boundaryBuffer[boundaryIndex[label]] = i; // 记录边界点，一个记录所有边界点
			boundaryIndex[label]++;// 记录索引，根据label可以索引到各个多边形的边界点
		}
	}

	//***********************************************************************
	//***********************************************************************

	//store the boundary locations stored by boundaryBuffer into
	//the region list for each region

	//***********************************************************************
	//***********************************************************************

	//destroy the old region list
	if(regionList)	delete regionList;

	//create a new region list
	int totalRegionPtCount = L;
	if(!(regionList	= new RegionList(regionCount, totalBoundaryCount, N, totalRegionPtCount)))
		ErrorHandler("msImageProcessor", "DefineRegions", "Not enough memory.");

	/************************************************************************/
	/* 记录每个多边形所包含的像素                                                                     */
	/************************************************************************/
	int	*regionBuffer	= new int [totalRegionPtCount];
	int *regionPtsIdx = new int[regionCount]; 
	regionPtsIdx[0]=0;
	for (j = 1; j < regionCount; j ++)
	{
		regionPtsIdx[j]=regionPtsIdx[j-1]+modePointCounts[j-1];
	}
	for (i = 0; i < L; i ++)
	{
		int label = labels[i];
		regionBuffer[regionPtsIdx[label]]=i;
		regionPtsIdx[label]++;
	}
	//add boundary locations for each region using the boundary
	//buffer and boundary counts
	counter	= 0;
	int counter2 = 0;
	for(i = 0; i < regionCount; i++)
	{
		/************************************************************************/
		/*    int label, int boundaryPtCount, int *boundaryIndeces ,int regionPtCount, int *regionIndeces                                                                  */
		/************************************************************************/
		regionList->AddRegion(i, boundaryCount[i], &boundaryBuffer[counter],modePointCounts[i],&regionBuffer[counter2]);
		counter += boundaryCount[i];
		counter2 += modePointCounts[i];
	}


	//***********************************************************************
	//***********************************************************************

   // dealocate local used memory
 	delete [] boundaryMap;// for boundary
   delete [] boundaryCount;// for RegionList 中region的boundarypointCount
	delete [] boundaryBuffer;// for RegionList 中indexTable
   delete [] boundaryIndex;
   delete [] regionPtsIdx;
   delete [] regionBuffer;

	//done.
	return;

}

	/*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
	/*  Image Data Searching/Distance Calculation */
	/*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/
		
/*******************************************************/
/*In Window                                            */
/*******************************************************/
/*Returns true if the two specified data points are    */
/*within rR of each other.                             */
/*******************************************************/
/*Pre:                                                 */
/*      - mode1 and mode2 are indeces into msRawData   */
/*        specifying the modes of the pixels having    */
/*        these indeces.                               */
/*Post:                                                */
/*      - true is returned if mode1 and mode2 are wi-  */
/*        thin rR of one another, false is returned    */
/*        otherwise.                                   */
/*******************************************************/

bool msImageProcessor::InWindow(int mode1, int mode2)
{
	int		k		= 1, s	= 0, p;
	double	diff	= 0, el;
	while((diff < 0.25)&&(k != kp)) // Partial Distortion Search
	{
		//Calculate distance squared of sub-space s	
		diff = 0;
		for(p = 0; p < P[k]; p++)
		{
			el    = (modes[mode1*N+p+s]-modes[mode2*N+p+s])/(h[k]*offset[k]);
			if((!p)&&(k == 1)&&(modes[mode1*N] > 80))
				diff += 4*el*el;
			else
				diff += el*el;
		}
		
		//next subspace
		s += P[k];
		k++;
	}
	return (bool)(diff < 0.25);
}

/*******************************************************/
/*Square Distance                                      */
/*******************************************************/
/*Computs the normalized square distance between two   */
/*modes.                                               */
/*******************************************************/
/*Pre:                                                 */
/*      - mode1 and mode2 are indeces into the modes   */
/*        array specifying two modes of the image      */
/*Post:                                                */
/*      - the normalized square distance between modes */
/*        indexed by mode1 and mode2 has been calc-    */
/*        ulated and the result has been returned.     */
/*******************************************************/

float msImageProcessor::SqDistance(int mode1, int mode2)
{

	int		k		= 1, s	= 0, p;
	float	dist	= 0, el;
	for(k = 1; k < kp; k++)
	{
		//Calculate distance squared of sub-space s	
		for(p = 0; p < P[k]; p++)
		{
			el    = (modes[mode1*N+p+s]-modes[mode2*N+p+s])/(h[k]*offset[k]);
			dist += el*el;
		}
		
		//next subspace
		s += P[k];
		k++;
	}

	//return normalized square distance between modes
	//1 and 2
	return dist;

}

	/*/\/\/\/\/\/\/\/\/\/\*/
	/*  Memory Management */
	/*\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Initialize Output                                    */
/*******************************************************/
/*Allocates memory needed by the mean shift image pro- */
/*cessor class output storage data structure.          */
/*******************************************************/
/*Post:                                                */
/*      - the memory needed by the output storage      */
/*        structure of this class has been (re-)allo-  */
/*        cated.                                       */
/*******************************************************/

void msImageProcessor::InitializeOutput( void )
{

	//De-allocate memory if output was defined for previous image
	DestroyOutput();

	//Allocate memory for msRawData (filtered image output)
	if(!(msRawData = new float [L*N]))
	{
		ErrorHandler("msImageProcessor", "Allocate", "Not enough memory.");
		return;
	}

	//Allocate memory used to store image modes and their corresponding regions...
	if((!(modes = new float [L*(N+2)]))||(!(labels = new int [L]))||(!(modePointCounts = new int [L]))||(!(indexTable	= new int [L])))
	{
		ErrorHandler("msImageProcessor", "Allocate", "Not enough memory");
		return;
	}

	//Allocate memory for integer modes used to perform connected components
	//(image labeling)...
//	if(!(LUV_data = new	int [N*L]))
   if (!(LUV_data = new float[N*L]))
	{
		ErrorHandler("msImageProcessor", "Allocate", "Not enough memory");
		return;
	}

	//indicate that the class output storage structure has been defined
	class_state.OUTPUT_DEFINED	= true;

}

/*******************************************************/
/*Destroy Output                                       */
/*******************************************************/
/*De-allocates memory needed by the mean shift image   */
/*processor class output storage data structure.       */
/*******************************************************/
/*Post:                                                */
/*      - the memory needed by the output storage      */
/*        structure of this class has been de-alloc-   */
/*        ated.                                        */
/*      - the output storage structure has been init-  */
/*        ialized for re-use.                          */
/*******************************************************/

void msImageProcessor::DestroyOutput( void )
{

	//de-allocate memory for msRawData (filtered image output)
	if (msRawData)			delete [] msRawData;

	//de-allocate memory used by output storage and image
	//classification structure
	if (modes)				delete [] modes;
	if (labels)				delete [] labels;
	if (modePointCounts)	delete [] modePointCounts;
	if (indexTable)			delete [] indexTable;
	
	//de-allocate memory for LUV_data
	if (LUV_data)			delete [] LUV_data;
		
	//initialize data members for re-use...

	//initialize output structures...
	msRawData			= NULL;

	//re-initialize classification structure
	modes						= NULL;
	labels						= NULL;
	modePointCounts				= NULL;
	regionCount					= 0;

	//indicate that the output has been destroyed
	class_state.OUTPUT_DEFINED	= false;

	//done.
	return;

}

// NEW
void msImageProcessor::NewOptimizedFilter1(float sigmaS, float sigmaR)
{
	// Declare Variables
	int		iterationCount, i, j, k, modeCandidateX, modeCandidateY, modeCandidate_i;
	double	mvAbs, diff, el;
	
	//make sure that a lattice height and width have
	//been defined...
	if(!height)
	{
		ErrorHandler("msImageProcessor", "LFilter", "Lattice height and width are undefined.");
		return;
	}

	//re-assign bandwidths to sigmaS and sigmaR
	if(((h[0] = sigmaS) <= 0)||((h[1] = sigmaR) <= 0))
	{
		ErrorHandler("msImageProcessor", "Segment", "sigmaS and/or sigmaR is zero or negative.");
		return;
	}
	
	//define input data dimension with lattice
	int lN	= N + 2;
	
	// Traverse each data point applying mean shift
	// to each data point
	
	// Allcocate memory for yk
	double	*yk		= new double [lN];
	
	// Allocate memory for Mh
	double	*Mh		= new double [lN];

   // let's use some temporary data
   float* sdata;
   sdata = new float[lN*L];

   // copy the scaled data
   int idxs, idxd;
   idxs = idxd = 0;
   if (N==3)
   {
      for(i=0; i<L; i++)
      {
         sdata[idxs++] = (i%width)/sigmaS;
         sdata[idxs++] = (i/width)/sigmaS;
         sdata[idxs++] = data[idxd++]/sigmaR;
         sdata[idxs++] = data[idxd++]/sigmaR;
         sdata[idxs++] = data[idxd++]/sigmaR;
      }
   } else if (N==1)
   {
      for(i=0; i<L; i++)
      {
         sdata[idxs++] = (i%width)/sigmaS;
         sdata[idxs++] = (i/width)/sigmaS;
         sdata[idxs++] = data[idxd++]/sigmaR;
      }
   } else
   {
      for(i=0; i<L; i++)
      {
         sdata[idxs++] = (i%width)/sigmaS;
         sdata[idxs++] = (i/width)/sigmaS;
         for (j=0; j<N; j++)
            sdata[idxs++] = data[idxd++]/sigmaR;
      }
   }
   // index the data in the 3d buckets (x, y, L)
   int* buckets;
   int* slist;
   slist = new int[L];
   int bucNeigh[27];

   float sMins; // just for L
   float sMaxs[3]; // for all
   sMaxs[0] = width/sigmaS;
   sMaxs[1] = height/sigmaS;
   sMins = sMaxs[2] = sdata[2];
   idxs = 2;
   float cval;
   for(i=0; i<L; i++)
   {
      cval = sdata[idxs];
      if (cval < sMins)
         sMins = cval;
      else if (cval > sMaxs[2])
         sMaxs[2] = cval;

      idxs += lN;
   }

   int nBuck1, nBuck2, nBuck3;
   int cBuck1, cBuck2, cBuck3, cBuck;
   nBuck1 = (int) (sMaxs[0] + 3);
   nBuck2 = (int) (sMaxs[1] + 3);
   nBuck3 = (int) (sMaxs[2] - sMins + 3);
   buckets = new int[nBuck1*nBuck2*nBuck3];
   for(i=0; i<(nBuck1*nBuck2*nBuck3); i++)
      buckets[i] = -1;

   idxs = 0;
   for(i=0; i<L; i++)
   {
      // find bucket for current data and add it to the list
      cBuck1 = (int) sdata[idxs] + 1;
      cBuck2 = (int) sdata[idxs+1] + 1;
      cBuck3 = (int) (sdata[idxs+2] - sMins) + 1;
      cBuck = cBuck1 + nBuck1*(cBuck2 + nBuck2*cBuck3);

      slist[i] = buckets[cBuck];
      buckets[cBuck] = i;

      idxs += lN;
   }
   // init bucNeigh
   idxd = 0;
   for (cBuck1=-1; cBuck1<=1; cBuck1++)
   {
      for (cBuck2=-1; cBuck2<=1; cBuck2++)
      {
         for (cBuck3=-1; cBuck3<=1; cBuck3++)
         {
            bucNeigh[idxd++] = cBuck1 + nBuck1*(cBuck2 + nBuck2*cBuck3);
         }
      }
   }
   double wsuml, weight;
   double hiLTr = 80.0/sigmaR;
   // done indexing/hashing

	
	// Initialize mode table used for basin of attraction
	memset(modeTable, 0, width*height);
	
	// proceed ...
#ifdef PROMPT
	msSys.Prompt("done.\nApplying mean shift (Using Lattice) ... ");
#ifdef SHOW_PROGRESS
	msSys.Prompt("\n 0%%");
#endif
#endif


	for(i = 0; i < L; i++)
	{
		// if a mode was already assigned to this data point
		// then skip this point, otherwise proceed to
		// find its mode by applying mean shift...
		if (modeTable[i] == 1)
			continue;

		// initialize point list...
		pointCount = 0;

		// Assign window center (window centers are
		// initialized by createLattice to be the point
		// data[i])
      idxs = i*lN;
      for (j=0; j<lN; j++)
         yk[j] = sdata[idxs+j];
		
		// Calculate the mean shift vector using the lattice
		// LatticeMSVector(Mh, yk); // modify to new
      /*****************************************************/
   	// Initialize mean shift vector
	   for(j = 0; j < lN; j++)
   		Mh[j] = 0;
   	wsuml = 0;
      // uniformLSearch(Mh, yk_ptr); // modify to new
      // find bucket of yk
      cBuck1 = (int) yk[0] + 1;
      cBuck2 = (int) yk[1] + 1;
      cBuck3 = (int) (yk[2] - sMins) + 1;
      cBuck = cBuck1 + nBuck1*(cBuck2 + nBuck2*cBuck3);
      for (j=0; j<27; j++)
      {
         idxd = buckets[cBuck+bucNeigh[j]];
         // list parse, crt point is cHeadList
         while (idxd>=0)
         {
            idxs = lN*idxd;
            // determine if inside search window
            el = sdata[idxs+0]-yk[0];
            diff = el*el;
            el = sdata[idxs+1]-yk[1];
            diff += el*el;

            if (diff < 1.0)
            {
               el = sdata[idxs+2]-yk[2];
               if (yk[2] > hiLTr)
                  diff = 4*el*el;
               else
                  diff = el*el;

               if (N>1)
               {
                  el = sdata[idxs+3]-yk[3];
                  diff += el*el;
                  el = sdata[idxs+4]-yk[4];
                  diff += el*el;
               }

               if (diff < 1.0)
               {
                  weight = 1-weightMap[idxd];
                  for (k=0; k<lN; k++)
                     Mh[k] += weight*sdata[idxs+k];
                  wsuml += weight;
               }
            }
            idxd = slist[idxd];
         }
      }
   	if (wsuml > 0)
   	{
		   for(j = 0; j < lN; j++)
   			Mh[j] = Mh[j]/wsuml - yk[j];
   	}
   	else
   	{
		   for(j = 0; j < lN; j++)
   			Mh[j] = 0;
   	}
      /*****************************************************/
   	// Calculate its magnitude squared
		//mvAbs = 0;
		//for(j = 0; j < lN; j++)
		//	mvAbs += Mh[j]*Mh[j];
      mvAbs = (Mh[0]*Mh[0]+Mh[1]*Mh[1])*sigmaS*sigmaS;
      if (N==3)
         mvAbs += (Mh[2]*Mh[2]+Mh[3]*Mh[3]+Mh[4]*Mh[4])*sigmaR*sigmaR;
      else
         mvAbs += Mh[2]*Mh[2]*sigmaR*sigmaR;

		
		// Keep shifting window center until the magnitude squared of the
		// mean shift vector calculated at the window center location is
		// under a specified threshold (Epsilon)
		
		// NOTE: iteration count is for speed up purposes only - it
		//       does not have any theoretical importance
		iterationCount = 1;
		while((mvAbs >= EPSILON)&&(iterationCount < LIMIT))
		{
			
			// Shift window location
			for(j = 0; j < lN; j++)
				yk[j] += Mh[j];
			
			// check to see if the current mode location is in the
			// basin of attraction...

			// calculate the location of yk on the lattice
			modeCandidateX	= (int) (sigmaS*yk[0]+0.5);
			modeCandidateY	= (int) (sigmaS*yk[1]+0.5);
			modeCandidate_i	= modeCandidateY*width + modeCandidateX;

			// if mvAbs != 0 (yk did indeed move) then check
			// location basin_i in the mode table to see if
			// this data point either:
			
			// (1) has not been associated with a mode yet
			//     (modeTable[basin_i] = 0), so associate
			//     it with this one
			//
			// (2) it has been associated with a mode other
			//     than the one that this data point is converging
			//     to (modeTable[basin_i] = 1), so assign to
			//     this data point the same mode as that of basin_i

			if ((modeTable[modeCandidate_i] != 2) && (modeCandidate_i != i))
			{
				// obtain the data point at basin_i to
				// see if it is within h*TC_DIST_FACTOR of
				// of yk
            diff = 0;
            idxs = lN*modeCandidate_i;
            for (k=2; k<lN; k++)
            {
               el = sdata[idxs+k] - yk[k];
               diff += el*el;
            }

				// if the data point at basin_i is within
				// a distance of h*TC_DIST_FACTOR of yk
				// then depending on modeTable[basin_i] perform
				// either (1) or (2)
				if (diff < TC_DIST_FACTOR)
				{
					// if the data point at basin_i has not
					// been associated to a mode then associate
					// it with the mode that this one will converge
					// to
					if (modeTable[modeCandidate_i] == 0)
					{
						// no mode associated yet so associate
						// it with this one...
						pointList[pointCount++]		= modeCandidate_i;
						modeTable[modeCandidate_i]	= 2;

					} else
					{

						// the mode has already been associated with
						// another mode, thererfore associate this one
						// mode and the modes in the point list with
						// the mode associated with data[basin_i]...

						// store the mode info into yk using msRawData...
						for (j = 0; j < N; j++)
							yk[j+2] = msRawData[modeCandidate_i*N+j]/sigmaR;

						// update mode table for this data point
						// indicating that a mode has been associated
						// with it
						modeTable[i] = 1;

						// indicate that a mode has been associated
						// to this data point (data[i])
						mvAbs = -1;

						// stop mean shift calculation...
						break;
					}
				}
			}
			
         // Calculate the mean shift vector at the new
         // window location using lattice
         // Calculate the mean shift vector using the lattice
         // LatticeMSVector(Mh, yk); // modify to new
         /*****************************************************/
         // Initialize mean shift vector
         for(j = 0; j < lN; j++)
            Mh[j] = 0;
         wsuml = 0;
         // uniformLSearch(Mh, yk_ptr); // modify to new
         // find bucket of yk
         cBuck1 = (int) yk[0] + 1;
         cBuck2 = (int) yk[1] + 1;
         cBuck3 = (int) (yk[2] - sMins) + 1;
         cBuck = cBuck1 + nBuck1*(cBuck2 + nBuck2*cBuck3);
         for (j=0; j<27; j++)
         {
            idxd = buckets[cBuck+bucNeigh[j]];
            // list parse, crt point is cHeadList
            while (idxd>=0)
            {
               idxs = lN*idxd;
               // determine if inside search window
               el = sdata[idxs+0]-yk[0];
               diff = el*el;
               el = sdata[idxs+1]-yk[1];
               diff += el*el;
               
               if (diff < 1.0)
               {
                  el = sdata[idxs+2]-yk[2];
                  if (yk[2] > hiLTr)
                     diff = 4*el*el;
                  else
                     diff = el*el;
                  
                  if (N>1)
                  {
                     el = sdata[idxs+3]-yk[3];
                     diff += el*el;
                     el = sdata[idxs+4]-yk[4];
                     diff += el*el;
                  }
                  
                  if (diff < 1.0)
                  {
                     weight = 1-weightMap[idxd];
                     for (k=0; k<lN; k++)
                        Mh[k] += weight*sdata[idxs+k];
                     wsuml += weight;
                  }
               }
               idxd = slist[idxd];
            }
         }
         if (wsuml > 0)
         {
            for(j = 0; j < lN; j++)
               Mh[j] = Mh[j]/wsuml - yk[j];
         }
         else
         {
            for(j = 0; j < lN; j++)
               Mh[j] = 0;
         }
         /*****************************************************/
			
			// Calculate its magnitude squared
			//mvAbs = 0;
			//for(j = 0; j < lN; j++)
			//	mvAbs += Mh[j]*Mh[j];
         mvAbs = (Mh[0]*Mh[0]+Mh[1]*Mh[1])*sigmaS*sigmaS;
         if (N==3)
            mvAbs += (Mh[2]*Mh[2]+Mh[3]*Mh[3]+Mh[4]*Mh[4])*sigmaR*sigmaR;
         else
            mvAbs += Mh[2]*Mh[2]*sigmaR*sigmaR;

			// Increment iteration count
			iterationCount++;
			
		}

		// if a mode was not associated with this data point
		// yet associate it with yk...
		if (mvAbs >= 0)
		{
			// Shift window location
			for(j = 0; j < lN; j++)
				yk[j] += Mh[j];
			
			// update mode table for this data point
			// indicating that a mode has been associated
			// with it
			modeTable[i] = 1;

		}
		
      for (k=0; k<N; k++)
         yk[k+2] *= sigmaR;

		// associate the data point indexed by
		// the point list with the mode stored
		// by yk
		for (j = 0; j < pointCount; j++)
		{
			// obtain the point location from the
			// point list
			modeCandidate_i = pointList[j];

			// update the mode table for this point
			modeTable[modeCandidate_i] = 1;

			//store result into msRawData...
			for(k = 0; k < N; k++)
				msRawData[N*modeCandidate_i+k] = (float)(yk[k+2]);
		}

		//store result into msRawData...
		for(j = 0; j < N; j++)
			msRawData[N*i+j] = (float)(yk[j+2]);

		// Prompt user on progress
#ifdef SHOW_PROGRESS
		percent_complete = (float)(i/(float)(L))*100;
		msSys.Prompt("\r%2d%%", (int)(percent_complete + 0.5));
#endif
	
		// Check to see if the algorithm has been halted
		if((i%PROGRESS_RATE == 0)&&((ErrorStatus = msSys.Progress((float)(i/(float)(L))*(float)(0.8)))) == EL_HALT)
			break;		
	}
	
	// Prompt user that filtering is completed
#ifdef PROMPT
#ifdef SHOW_PROGRESS
	msSys.Prompt("\r");
#endif
	msSys.Prompt("done.");
#endif
	// de-allocate memory
   delete [] buckets;
   delete [] slist;
   delete [] sdata;

	delete [] yk;
	delete [] Mh;
	
	// done.
	return;

}

// NEW
void msImageProcessor::NewOptimizedFilter2(float sigmaS, float sigmaR)
{
	// Declare Variables
	int		iterationCount, i, j, k, modeCandidateX, modeCandidateY, modeCandidate_i;
	double	mvAbs, diff, el;
	
	//make sure that a lattice height and width have
	//been defined...
	if(!height)
	{
		ErrorHandler("msImageProcessor", "LFilter", "Lattice height and width are undefined.");
		return;
	}

	//re-assign bandwidths to sigmaS and sigmaR
	if(((h[0] = sigmaS) <= 0)||((h[1] = sigmaR) <= 0))
	{
		ErrorHandler("msImageProcessor", "Segment", "sigmaS and/or sigmaR is zero or negative.");
		return;
	}
	
	//define input data dimension with lattice
	int lN	= N + 2;
	
	// Traverse each data point applying mean shift
	// to each data point
	
	// Allcocate memory for yk
	double	*yk		= new double [lN];
	
	// Allocate memory for Mh
	double	*Mh		= new double [lN];

   // let's use some temporary data
   float* sdata;
   sdata = new float[lN*L];

   // copy the scaled data
   int idxs, idxd;
   idxs = idxd = 0;
   if (N==3)
   {
      for(i=0; i<L; i++)
      {
         sdata[idxs++] = (i%width)/sigmaS;
         sdata[idxs++] = (i/width)/sigmaS;
         sdata[idxs++] = data[idxd++]/sigmaR;
         sdata[idxs++] = data[idxd++]/sigmaR;
         sdata[idxs++] = data[idxd++]/sigmaR;
      }
   } else if (N==1)
   {
      for(i=0; i<L; i++)
      {
         sdata[idxs++] = (i%width)/sigmaS;
         sdata[idxs++] = (i/width)/sigmaS;
         sdata[idxs++] = data[idxd++]/sigmaR;
      }
   } else
   {
      for(i=0; i<L; i++)
      {
         sdata[idxs++] = (i%width)/sigmaS;
         sdata[idxs++] = (i/width)/sigmaS;
         for (j=0; j<N; j++)
            sdata[idxs++] = data[idxd++]/sigmaR;
      }
   }
   // index the data in the 3d buckets (x, y, L)
   int* buckets;
   int* slist;
   slist = new int[L];
   int bucNeigh[27];

   float sMins; // just for L
   float sMaxs[3]; // for all
   sMaxs[0] = width/sigmaS;
   sMaxs[1] = height/sigmaS;
   sMins = sMaxs[2] = sdata[2];
   idxs = 2;
   float cval;
   for(i=0; i<L; i++)
   {
      cval = sdata[idxs];
      if (cval < sMins)
         sMins = cval;
      else if (cval > sMaxs[2])
         sMaxs[2] = cval;

      idxs += lN;
   }

   int nBuck1, nBuck2, nBuck3;
   int cBuck1, cBuck2, cBuck3, cBuck;
   nBuck1 = (int) (sMaxs[0] + 3);
   nBuck2 = (int) (sMaxs[1] + 3);
   nBuck3 = (int) (sMaxs[2] - sMins + 3);
   buckets = new int[nBuck1*nBuck2*nBuck3];
   for(i=0; i<(nBuck1*nBuck2*nBuck3); i++)
      buckets[i] = -1;

   idxs = 0;
   for(i=0; i<L; i++)
   {
      // find bucket for current data and add it to the list
      cBuck1 = (int) sdata[idxs] + 1;
      cBuck2 = (int) sdata[idxs+1] + 1;
      cBuck3 = (int) (sdata[idxs+2] - sMins) + 1;
      cBuck = cBuck1 + nBuck1*(cBuck2 + nBuck2*cBuck3);

      slist[i] = buckets[cBuck];
      buckets[cBuck] = i;

      idxs += lN;
   }
   // init bucNeigh
   idxd = 0;
   for (cBuck1=-1; cBuck1<=1; cBuck1++)
   {
      for (cBuck2=-1; cBuck2<=1; cBuck2++)
      {
         for (cBuck3=-1; cBuck3<=1; cBuck3++)
         {
            bucNeigh[idxd++] = cBuck1 + nBuck1*(cBuck2 + nBuck2*cBuck3);
         }
      }
   }
   double wsuml, weight;
   double hiLTr = 80.0/sigmaR;
   // done indexing/hashing

	
	// Initialize mode table used for basin of attraction
	memset(modeTable, 0, width*height);
	
	// proceed ...
#ifdef PROMPT
	msSys.Prompt("done.\nApplying mean shift (Using Lattice) ... ");
#ifdef SHOW_PROGRESS
	msSys.Prompt("\n 0%%");
#endif
#endif


	for(i = 0; i < L; i++)
	{
		// if a mode was already assigned to this data point
		// then skip this point, otherwise proceed to
		// find its mode by applying mean shift...
		if (modeTable[i] == 1)
			continue;

		// initialize point list...
		pointCount = 0;

		// Assign window center (window centers are
		// initialized by createLattice to be the point
		// data[i])
      idxs = i*lN;
      for (j=0; j<lN; j++)
         yk[j] = sdata[idxs+j];
		
		// Calculate the mean shift vector using the lattice
		// LatticeMSVector(Mh, yk); // modify to new
      /*****************************************************/
   	// Initialize mean shift vector
	   for(j = 0; j < lN; j++)
   		Mh[j] = 0;
   	wsuml = 0;
      // uniformLSearch(Mh, yk_ptr); // modify to new
      // find bucket of yk
      cBuck1 = (int) yk[0] + 1;
      cBuck2 = (int) yk[1] + 1;
      cBuck3 = (int) (yk[2] - sMins) + 1;
      cBuck = cBuck1 + nBuck1*(cBuck2 + nBuck2*cBuck3);
      for (j=0; j<27; j++)
      {
         idxd = buckets[cBuck+bucNeigh[j]];
         // list parse, crt point is cHeadList
         while (idxd>=0)
         {
            idxs = lN*idxd;
            // determine if inside search window
            el = sdata[idxs+0]-yk[0];
            diff = el*el;
            el = sdata[idxs+1]-yk[1];
            diff += el*el;

            if (diff < 1.0)
            {
               el = sdata[idxs+2]-yk[2];
               if (yk[2] > hiLTr)
                  diff = 4*el*el;
               else
                  diff = el*el;

               if (N>1)
               {
                  el = sdata[idxs+3]-yk[3];
                  diff += el*el;
                  el = sdata[idxs+4]-yk[4];
                  diff += el*el;
               }

               if (diff < 1.0)
               {
                  weight = 1-weightMap[idxd];
                  for (k=0; k<lN; k++)
                     Mh[k] += weight*sdata[idxs+k];
                  wsuml += weight;

      				//set basin of attraction mode table
                  if (diff < speedThreshold)
                  {
				         if(modeTable[idxd] == 0)
				         {
         					pointList[pointCount++]	= idxd;
					         modeTable[idxd]	= 2;
      				   }
                  }
               }
            }
            idxd = slist[idxd];
         }
      }
   	if (wsuml > 0)
   	{
		   for(j = 0; j < lN; j++)
   			Mh[j] = Mh[j]/wsuml - yk[j];
   	}
   	else
   	{
		   for(j = 0; j < lN; j++)
   			Mh[j] = 0;
   	}
      /*****************************************************/
   	// Calculate its magnitude squared
		//mvAbs = 0;
		//for(j = 0; j < lN; j++)
		//	mvAbs += Mh[j]*Mh[j];
      mvAbs = (Mh[0]*Mh[0]+Mh[1]*Mh[1])*sigmaS*sigmaS;
      if (N==3)
         mvAbs += (Mh[2]*Mh[2]+Mh[3]*Mh[3]+Mh[4]*Mh[4])*sigmaR*sigmaR;
      else
         mvAbs += Mh[2]*Mh[2]*sigmaR*sigmaR;

		
		// Keep shifting window center until the magnitude squared of the
		// mean shift vector calculated at the window center location is
		// under a specified threshold (Epsilon)
		
		// NOTE: iteration count is for speed up purposes only - it
		//       does not have any theoretical importance
		iterationCount = 1;
		while((mvAbs >= EPSILON)&&(iterationCount < LIMIT))
		{
			
			// Shift window location
			for(j = 0; j < lN; j++)
				yk[j] += Mh[j];
			
			// check to see if the current mode location is in the
			// basin of attraction...

			// calculate the location of yk on the lattice
			modeCandidateX	= (int) (sigmaS*yk[0]+0.5);
			modeCandidateY	= (int) (sigmaS*yk[1]+0.5);
			modeCandidate_i	= modeCandidateY*width + modeCandidateX;

			// if mvAbs != 0 (yk did indeed move) then check
			// location basin_i in the mode table to see if
			// this data point either:
			
			// (1) has not been associated with a mode yet
			//     (modeTable[basin_i] = 0), so associate
			//     it with this one
			//
			// (2) it has been associated with a mode other
			//     than the one that this data point is converging
			//     to (modeTable[basin_i] = 1), so assign to
			//     this data point the same mode as that of basin_i

			if ((modeTable[modeCandidate_i] != 2) && (modeCandidate_i != i))
			{
				// obtain the data point at basin_i to
				// see if it is within h*TC_DIST_FACTOR of
				// of yk
            diff = 0;
            idxs = lN*modeCandidate_i;
            for (k=2; k<lN; k++)
            {
               el = sdata[idxs+k] - yk[k];
               diff += el*el;
            }

				// if the data point at basin_i is within
				// a distance of h*TC_DIST_FACTOR of yk
				// then depending on modeTable[basin_i] perform
				// either (1) or (2)
				if (diff < speedThreshold)
				{
					// if the data point at basin_i has not
					// been associated to a mode then associate
					// it with the mode that this one will converge
					// to
					if (modeTable[modeCandidate_i] == 0)
					{
						// no mode associated yet so associate
						// it with this one...
						pointList[pointCount++]		= modeCandidate_i;
						modeTable[modeCandidate_i]	= 2;

					} else
					{

						// the mode has already been associated with
						// another mode, thererfore associate this one
						// mode and the modes in the point list with
						// the mode associated with data[basin_i]...

						// store the mode info into yk using msRawData...
						for (j = 0; j < N; j++)
							yk[j+2] = msRawData[modeCandidate_i*N+j]/sigmaR;

						// update mode table for this data point
						// indicating that a mode has been associated
						// with it
						modeTable[i] = 1;

						// indicate that a mode has been associated
						// to this data point (data[i])
						mvAbs = -1;

						// stop mean shift calculation...
						break;
					}
				}
			}
			
         // Calculate the mean shift vector at the new
         // window location using lattice
         // Calculate the mean shift vector using the lattice
         // LatticeMSVector(Mh, yk); // modify to new
         /*****************************************************/
         // Initialize mean shift vector
         for(j = 0; j < lN; j++)
            Mh[j] = 0;
         wsuml = 0;
         // uniformLSearch(Mh, yk_ptr); // modify to new
         // find bucket of yk
         cBuck1 = (int) yk[0] + 1;
         cBuck2 = (int) yk[1] + 1;
         cBuck3 = (int) (yk[2] - sMins) + 1;
         cBuck = cBuck1 + nBuck1*(cBuck2 + nBuck2*cBuck3);
         for (j=0; j<27; j++)
         {
            idxd = buckets[cBuck+bucNeigh[j]];
            // list parse, crt point is cHeadList
            while (idxd>=0)
            {
               idxs = lN*idxd;
               // determine if inside search window
               el = sdata[idxs+0]-yk[0];
               diff = el*el;
               el = sdata[idxs+1]-yk[1];
               diff += el*el;
               
               if (diff < 1.0)
               {
                  el = sdata[idxs+2]-yk[2];
                  if (yk[2] > hiLTr)
                     diff = 4*el*el;
                  else
                     diff = el*el;
                  
                  if (N>1)
                  {
                     el = sdata[idxs+3]-yk[3];
                     diff += el*el;
                     el = sdata[idxs+4]-yk[4];
                     diff += el*el;
                  }
                  
                  if (diff < 1.0)
                  {
                     weight = 1-weightMap[idxd];
                     for (k=0; k<lN; k++)
                        Mh[k] += weight*sdata[idxs+k];
                     wsuml += weight;

         				//set basin of attraction mode table
                     if (diff < speedThreshold)
                     {
   				         if(modeTable[idxd] == 0)
				            {
            					pointList[pointCount++]	= idxd;
					            modeTable[idxd]	= 2;
      				      }
                     }

                  }
               }
               idxd = slist[idxd];
            }
         }
         if (wsuml > 0)
         {
            for(j = 0; j < lN; j++)
               Mh[j] = Mh[j]/wsuml - yk[j];
         }
         else
         {
            for(j = 0; j < lN; j++)
               Mh[j] = 0;
         }
         /*****************************************************/
			
			// Calculate its magnitude squared
			//mvAbs = 0;
			//for(j = 0; j < lN; j++)
			//	mvAbs += Mh[j]*Mh[j];
         mvAbs = (Mh[0]*Mh[0]+Mh[1]*Mh[1])*sigmaS*sigmaS;
         if (N==3)
            mvAbs += (Mh[2]*Mh[2]+Mh[3]*Mh[3]+Mh[4]*Mh[4])*sigmaR*sigmaR;
         else
            mvAbs += Mh[2]*Mh[2]*sigmaR*sigmaR;

			// Increment iteration count
			iterationCount++;
			
		}

		// if a mode was not associated with this data point
		// yet associate it with yk...
		if (mvAbs >= 0)
		{
			// Shift window location
			for(j = 0; j < lN; j++)
				yk[j] += Mh[j];
			
			// update mode table for this data point
			// indicating that a mode has been associated
			// with it
			modeTable[i] = 1;

		}
		
      for (k=0; k<N; k++)
         yk[k+2] *= sigmaR;

		// associate the data point indexed by
		// the point list with the mode stored
		// by yk
		for (j = 0; j < pointCount; j++)
		{
			// obtain the point location from the
			// point list
			modeCandidate_i = pointList[j];

			// update the mode table for this point
			modeTable[modeCandidate_i] = 1;

			//store result into msRawData...
			for(k = 0; k < N; k++)
				msRawData[N*modeCandidate_i+k] = (float)(yk[k+2]);
		}

		//store result into msRawData...
		for(j = 0; j < N; j++)
			msRawData[N*i+j] = (float)(yk[j+2]);

		// Prompt user on progress
#ifdef SHOW_PROGRESS
		percent_complete = (float)(i/(float)(L))*100;
		msSys.Prompt("\r%2d%%", (int)(percent_complete + 0.5));
#endif
	
		// Check to see if the algorithm has been halted
		if((i%PROGRESS_RATE == 0)&&((ErrorStatus = msSys.Progress((float)(i/(float)(L))*(float)(0.8)))) == EL_HALT)
			break;		
	}
	
	// Prompt user that filtering is completed
#ifdef PROMPT
#ifdef SHOW_PROGRESS
	msSys.Prompt("\r");
#endif
	msSys.Prompt("done.");
#endif
	// de-allocate memory
   delete [] buckets;
   delete [] slist;
   delete [] sdata;

	delete [] yk;
	delete [] Mh;
	
	// done.
	return;

}

void msImageProcessor::NewNonOptimizedFilter(float sigmaS, float sigmaR)
{

	// Declare Variables
	int   iterationCount, i, j, k;
	double mvAbs, diff, el;
	
	//make sure that a lattice height and width have
	//been defined...
	if(!height)
	{
		ErrorHandler("msImageProcessor", "LFilter", "Lattice height and width are undefined.");
		return;
	}

	//re-assign bandwidths to sigmaS and sigmaR
	if(((h[0] = sigmaS) <= 0)||((h[1] = sigmaR) <= 0))
	{
		ErrorHandler("msImageProcessor", "Segment", "sigmaS and/or sigmaR is zero or negative.");
		return;
	}
	
	//define input data dimension with lattice
	int lN	= N + 2;
	
	// Traverse each data point applying mean shift
	// to each data point
	
	// Allcocate memory for yk
	double	*yk		= new double [lN];
	
	// Allocate memory for Mh
	double	*Mh		= new double [lN];

   // let's use some temporary data
   double* sdata;
   sdata = new double[lN*L];

   // copy the scaled data
   int idxs, idxd;
   idxs = idxd = 0;
   if (N==3)
   {
      for(i=0; i<L; i++)
      {
         sdata[idxs++] = (i%width)/sigmaS;
         sdata[idxs++] = (i/width)/sigmaS;
         sdata[idxs++] = data[idxd++]/sigmaR;
         sdata[idxs++] = data[idxd++]/sigmaR;
         sdata[idxs++] = data[idxd++]/sigmaR;
      }
   } 
   else if (N==1)
   {
      for(i=0; i<L; i++)
      {
         sdata[idxs++] = (i%width)/sigmaS;
         sdata[idxs++] = (i/width)/sigmaS;
         sdata[idxs++] = data[idxd++]/sigmaR;
      }
   }
   else
   {
      for(i=0; i<L; i++)
      {
         sdata[idxs++] = (i%width)/sigmaS;
         sdata[idxs++] = (i%width)/sigmaS;
         for (j=0; j<N; j++)
            sdata[idxs++] = data[idxd++]/sigmaR;
      }
   }
   // index the data in the 3d buckets (x, y, L)
   int* buckets;
   int* slist;
   slist = new int[L];
   int bucNeigh[27];

   double sMins; // just for L
   double sMaxs[3]; // for all
   sMaxs[0] = width/sigmaS;
   sMaxs[1] = height/sigmaS;
   sMins = sMaxs[2] = sdata[2];
   idxs = 2;
   double cval;
   for(i=0; i<L; i++)
   {
      cval = sdata[idxs];
      if (cval < sMins)
         sMins = cval;
      else if (cval > sMaxs[2])
         sMaxs[2] = cval;

      idxs += lN;
   }

   int nBuck1, nBuck2, nBuck3;
   int cBuck1, cBuck2, cBuck3, cBuck;
   nBuck1 = (int) (sMaxs[0] + 3);
   nBuck2 = (int) (sMaxs[1] + 3);
   nBuck3 = (int) (sMaxs[2] - sMins + 3);
   buckets = new int[nBuck1*nBuck2*nBuck3];
   for(i=0; i<(nBuck1*nBuck2*nBuck3); i++)
      buckets[i] = -1;

   idxs = 0;
   for(i=0; i<L; i++)
   {
      // find bucket for current data and add it to the list
      cBuck1 = (int) sdata[idxs] + 1;
      cBuck2 = (int) sdata[idxs+1] + 1;
      cBuck3 = (int) (sdata[idxs+2] - sMins) + 1;
      cBuck = cBuck1 + nBuck1*(cBuck2 + nBuck2*cBuck3);

      slist[i] = buckets[cBuck];
      buckets[cBuck] = i;

      idxs += lN;
   }
   // init bucNeigh
   idxd = 0;
   for (cBuck1=-1; cBuck1<=1; cBuck1++)
   {
      for (cBuck2=-1; cBuck2<=1; cBuck2++)
      {
         for (cBuck3=-1; cBuck3<=1; cBuck3++)
         {
            bucNeigh[idxd++] = cBuck1 + nBuck1*(cBuck2 + nBuck2*cBuck3);
         }
      }
   }
   double wsuml, weight;
   double hiLTr = 80.0/sigmaR;
   // done indexing/hashing
	
	// proceed ...
#ifdef PROMPT
	msSys.Prompt("done.\nApplying mean shift (Using Lattice)... ");
#ifdef SHOW_PROGRESS
	msSys.Prompt("\n 0%%");
#endif
#endif

	for(i = 0; i < L; i++)
	{

		// Assign window center (window centers are
		// initialized by createLattice to be the point
		// data[i])
      idxs = i*lN;
      for (j=0; j<lN; j++)
         yk[j] = sdata[idxs+j];
		
		// Calculate the mean shift vector using the lattice
		// LatticeMSVector(Mh, yk);
      /*****************************************************/
   	// Initialize mean shift vector
	   for(j = 0; j < lN; j++)
   		Mh[j] = 0;
   	wsuml = 0;
      // uniformLSearch(Mh, yk_ptr); // modify to new
      // find bucket of yk
      cBuck1 = (int) yk[0] + 1;
      cBuck2 = (int) yk[1] + 1;
      cBuck3 = (int) (yk[2] - sMins) + 1;
      cBuck = cBuck1 + nBuck1*(cBuck2 + nBuck2*cBuck3);
      for (j=0; j<27; j++)
      {
         idxd = buckets[cBuck+bucNeigh[j]];
         // list parse, crt point is cHeadList
         while (idxd>=0)
         {
            idxs = lN*idxd;
            // determine if inside search window
            el = sdata[idxs+0]-yk[0];
            diff = el*el;
            el = sdata[idxs+1]-yk[1];
            diff += el*el;

            if (diff < 1.0)
            {
               el = sdata[idxs+2]-yk[2];
               if (yk[2] > hiLTr)
                  diff = 4*el*el;
               else
                  diff = el*el;

               if (N>1)
               {
                  el = sdata[idxs+3]-yk[3];
                  diff += el*el;
                  el = sdata[idxs+4]-yk[4];
                  diff += el*el;
               }

               if (diff < 1.0)
               {
                  weight = 1-weightMap[idxd];
                  for (k=0; k<lN; k++)
                     Mh[k] += weight*sdata[idxs+k];
                  wsuml += weight;
               }
            }
            idxd = slist[idxd];
         }
      }
   	if (wsuml > 0)
   	{
		   for(j = 0; j < lN; j++)
   			Mh[j] = Mh[j]/wsuml - yk[j];
   	}
   	else
   	{
		   for(j = 0; j < lN; j++)
   			Mh[j] = 0;
   	}
      /*****************************************************/
		
		// Calculate its magnitude squared
		mvAbs = 0;
		for(j = 0; j < lN; j++)
			mvAbs += Mh[j]*Mh[j];
		
		// Keep shifting window center until the magnitude squared of the
		// mean shift vector calculated at the window center location is
		// under a specified threshold (Epsilon)
		
		// NOTE: iteration count is for speed up purposes only - it
		//       does not have any theoretical importance
		iterationCount = 1;
		while((mvAbs >= EPSILON)&&(iterationCount < LIMIT))
		{
			
			// Shift window location
			for(j = 0; j < lN; j++)
				yk[j] += Mh[j];
			
			// Calculate the mean shift vector at the new
			// window location using lattice
			// LatticeMSVector(Mh, yk);
         /*****************************************************/
         // Initialize mean shift vector
         for(j = 0; j < lN; j++)
            Mh[j] = 0;
         wsuml = 0;
         // uniformLSearch(Mh, yk_ptr); // modify to new
         // find bucket of yk
         cBuck1 = (int) yk[0] + 1;
         cBuck2 = (int) yk[1] + 1;
         cBuck3 = (int) (yk[2] - sMins) + 1;
         cBuck = cBuck1 + nBuck1*(cBuck2 + nBuck2*cBuck3);
         for (j=0; j<27; j++)
         {
            idxd = buckets[cBuck+bucNeigh[j]];
            // list parse, crt point is cHeadList
            while (idxd>=0)
            {
               idxs = lN*idxd;
               // determine if inside search window
               el = sdata[idxs+0]-yk[0];
               diff = el*el;
               el = sdata[idxs+1]-yk[1];
               diff += el*el;
               
               if (diff < 1.0)
               {
                  el = sdata[idxs+2]-yk[2];
                  if (yk[2] > hiLTr)
                     diff = 4*el*el;
                  else
                     diff = el*el;
                  
                  if (N>1)
                  {
                     el = sdata[idxs+3]-yk[3];
                     diff += el*el;
                     el = sdata[idxs+4]-yk[4];
                     diff += el*el;
                  }
                  
                  if (diff < 1.0)
                  {
                     weight = 1-weightMap[idxd];
                     for (k=0; k<lN; k++)
                        Mh[k] += weight*sdata[idxs+k];
                     wsuml += weight;
                  }
               }
               idxd = slist[idxd];
            }
         }
         if (wsuml > 0)
         {
            for(j = 0; j < lN; j++)
               Mh[j] = Mh[j]/wsuml - yk[j];
         }
         else
         {
            for(j = 0; j < lN; j++)
               Mh[j] = 0;
         }
         /*****************************************************/
			
			// Calculate its magnitude squared
			//mvAbs = 0;
			//for(j = 0; j < lN; j++)
			//	mvAbs += Mh[j]*Mh[j];
         mvAbs = (Mh[0]*Mh[0]+Mh[1]*Mh[1])*sigmaS*sigmaS;
         if (N==3)
            mvAbs += (Mh[2]*Mh[2]+Mh[3]*Mh[3]+Mh[4]*Mh[4])*sigmaR*sigmaR;
         else
            mvAbs += Mh[2]*Mh[2]*sigmaR*sigmaR;

			// Increment interation count
			iterationCount++;
		}

		// Shift window location
		for(j = 0; j < lN; j++)
			yk[j] += Mh[j];
		
		//store result into msRawData...
		for(j = 0; j < N; j++)
			msRawData[N*i+j] = (float)(yk[j+2]*sigmaR);

		// Prompt user on progress
#ifdef SHOW_PROGRESS
		percent_complete = (float)(i/(float)(L))*100;
		msSys.Prompt("\r%2d%%", (int)(percent_complete + 0.5));
#endif
	
		// Check to see if the algorithm has been halted
		if((i%PROGRESS_RATE == 0)&&((ErrorStatus = msSys.Progress((float)(i/(float)(L))*(float)(0.8)))) == EL_HALT)
			break;
	}
	
	// Prompt user that filtering is completed
#ifdef PROMPT
#ifdef SHOW_PROGRESS
	msSys.Prompt("\r");
#endif
	msSys.Prompt("done.");
#endif
	
	// de-allocate memory
   delete [] buckets;
   delete [] slist;
   delete [] sdata;

	delete [] yk;
	delete [] Mh;

	// done.
	return;

}

void msImageProcessor::SetSpeedThreshold(float speedUpThreshold)
{
   speedThreshold = speedUpThreshold;
}

int * msImageProcessor::GetLabels()
{
	if(class_state.OUTPUT_DEFINED == false)
		return NULL;
	int *labels_;
	if(!(labels_ = new int [L]))
	{
		ErrorHandler("msImageProcessor", "GetLabels", "Not enough memory.");
		return NULL;
	}
	for(int i = 0; i < L; i++)
		labels_[i] = labels[i];

	return labels_;
}

int *msImageProcessor::GetMPC()
{
	if(class_state.OUTPUT_DEFINED == false)
		return NULL;
	int *MPC_out_;
	if(!(MPC_out_ = new int [regionCount]))
	{
		ErrorHandler("msImageProcessor", "GetMPC", "Not enough memory.");
		return NULL;
	}

	for(int i = 0; i < regionCount; i++)
		MPC_out_[i]	= modePointCounts[i]; //存储每种模式中所包含的点数，一个多边形一个模态，存储多边形的面积

	return MPC_out_;
}

// 最好还返回数据的维度比较好
float *msImageProcessor::GetModeData()
{
	if(class_state.OUTPUT_DEFINED == false)
		return NULL;

	float *modes_;

	if(!(modes_	= new float [regionCount*N]))
	{
		ErrorHandler("msImageProcessor", "GetModeData", "Not enough memory.");
		return NULL;
	}

	// 返回多边形的特征数据
	for(int i = 0; i < regionCount*N; i++)
		modes_[i]	= modes[i]; 
	return modes_;
}

// 获取每个多边形的灰度值
byte * msImageProcessor::ModeData2RGB()
{
	byte *regionsRGB;
	float *modes_;
	modes_ = GetModeData();
	if (modes_)
	{
		if (!(regionsRGB = new BYTE [regionCount*N]))
		{
			ErrorHandler("msImageProcessor", "ModeData2RGB", "Not enough memory.");
			return NULL;
		}
		LUVtoRGB(modes_,regionsRGB);
		return regionsRGB;
	}
	else
	{
		return NULL;
	}
}

// 用一个数组记录每个多边形所包含的像素的索引
int* msImageProcessor::GetRegionPtIdxs()
{
// 	for (int i = 0; i < L; i ++)
// 	{
// 
// 	}
	return NULL;	
}