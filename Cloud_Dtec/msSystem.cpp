// msSystem.cpp: implementation of the msSystem class.
//
//////////////////////////////////////////////////////////////////////


#include "msSystem.h"
//include needed system libraries
#include	<time.h>
#include	<stdio.h>
#include	<stdarg.h>
#include	<stdlib.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

/*******************************************************/
/*Class Constructor                                    */
/*******************************************************/
/*Constructs a mean shift system object.               */
/*******************************************************/
/*Post:                                                */
/*      - an msSystem object has been properly init-   */
/*        ialized.                                     */
/*******************************************************/

msSystem::msSystem( void )
{

	//initialize currentTime
	currentTime = clock();

	//done.

}

/*******************************************************/
/*Class Destructor                                     */
/*******************************************************/
/*Destroys a mean shift system object.                 */
/*******************************************************/
/*Post:                                                */
/*      - an msSystem object has been properly dest-   */
/*        royed.                                       */
/*******************************************************/

msSystem::~msSystem( void )
{
	/* do nothing */
}

 /*/\/\/\/\/\/\/\/\/\*/
 /*** System Timer ***/
 /*\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Start Timer                                          */
/*******************************************************/
/*Sets the mean shift system time to the current       */
/*system time.                                         */
/*******************************************************/
/*Post:                                                */
/*      - the mean shift system time has been set to   */
/*        the current system time.                     */
/*******************************************************/

void msSystem::StartTimer( void )
{

	//set msSystem time to system time
	currentTime = clock();

	//done.
	return;

}

/*******************************************************/
/*Elapsed Time                                         */
/*******************************************************/
/*Returns the amount of time in seconds since the      */
/*mean shift system time was last set.                 */
/*******************************************************/
/*Post:                                                */
/*      - the amount of time in seconds since the mean */
/*        shift system time was last set is returned.  */
/*******************************************************/

double msSystem::ElapsedTime( void )
{

	//return the amount of time elapsed in seconds
	//since the msSystem time was last set...
	return ((double) (clock() - currentTime))/(CLOCKS_PER_SEC);

}

 /*/\/\/\/\/\/\/\/\/\/\*/
 /***  System Output ***/
 /*\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Prompt                                               */
/*******************************************************/
/*Output a text message to the user.                   */
/*******************************************************/
/*Pre:                                                 */
/*      - PromptStr is a string containing delimeters  */
/*        that is to be output to the user.            */
/*      - a variable set of arguments is also passed   */
/*        to this method that are used to replace      */
/*        the delimeters contained by PromptStr        */
/*Post:                                                */
/*      - the delimeters of PromptStr have been        */
/*        replaced accordingly using the variable      */
/*        set of arguments and the resulting string    */
/*        has been output to the user.                 */
/*******************************************************/

//extern void bgLogVar(const char *, va_list);

void msSystem::Prompt(const char *PromptStr, ...)
{

	//obtain argument list using ANSI standard...
	va_list	argList;
	va_start(argList, PromptStr);

	//print the output string to stderr using
	//vfprintf
	//bgLogVar(PromptStr, argList);
	va_end(argList);

	//done.
	return;

}

/*******************************************************/
/*Progress                                             */
/*******************************************************/
/*The progress of a specific algorithm of the mean     */
/*shift library is output to the user.                 */
/*******************************************************/
/*Pre:                                                 */
/*		- percentComplete indicates the percentage     */
/*		  of the algorithm that has executed and is    */
/*		  a floating point number from zero to one     */
/*Post:                                                */
/*      - the percentComplete has been noted by the    */
/*        interface (possibly to update a progress     */
/*        bar).                                        */
/*      - if the thread executing the mean shift code  */
/*        must halt execution msSYS_HALT is returned,  */
/*        msSYS_OK is returned otherwise               */
/*******************************************************/

///////////////////////////////////////////////////////////////////
//NOTE: This implementation is specific to EDISON. In order
//      for one to port the mean shift class to another project
//      or program one must re-implement this method.
///////////////////////////////////////////////////////////////////

//is set by the GUI when the user presses the Cancel button
//on the wxWindows progress modal window; this flag indicates
//to the mean shift library that it is to halt execution.
//This parameter is used and checked in the method
//BgMdiSegmentChild::OnSegment.

/////////////////////////**********此处修改过**********************************//
//extern bool stop_flag;

//is updated in the msSystem::Progress method and indicates to
//the wxWindows progress modal window the percent complete, such
//that it may update its progress bar accordingly; This parameter
//is used and checked in the method BgMdiSegmentChild::OnSegment.

/////////////////////////**********此处修改过**********************************//
//extern int	percentDone;

ErrorLevel msSystem::Progress(float percentComplete)
{
//	percentDone	= (int)(percentComplete*100);

	//check stop flag and return appropriate system state
	ErrorLevel		myState = EL_OKAY;
//	if(stop_flag)	myState	= EL_HALT;

	//done.
	return myState;

}