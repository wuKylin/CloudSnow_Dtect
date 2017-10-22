// msSystem.h: interface for the msSystem class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MSSYSTEM_H__D6BC486B_CCB1_4B5F_A3C6_903AE5D14518__INCLUDED_)
#define AFX_MSSYSTEM_H__D6BC486B_CCB1_4B5F_A3C6_903AE5D14518__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//Include standard mean shift library type definitions
#include	"tdef.h"

//Include standard libraries needed for msSystem prototype
#include	<time.h>

//extern void bgLogFile(const char*, ...);

//Mean Shify System class prototype
class msSystem {

 public:

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /* Class Constructor and Destructor */
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

  msSystem( void ); //Default Constructor
 ~msSystem( void ); //Class Destructor

 /*/\/\/\/\/\/\/\*/
 /* System Timer */
 /*\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|				  *  Start Timer  *                  |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Initializes the system timer. The timer object   |//
  //|   synthesized by this class is initialized during  |//
  //|   construction of the msSystem class to be the     |//
  //|   current time during construction.                |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

 void StartTimer( void );

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|				 *  Elapsed Time  *                  |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Returns the amount of time elapsed in seconds    |//
  //|   from when StartTimer() was called. If            |//
  //|   StartTimer() was not called, the time returned   |//
  //|   is the time elapsed from the construction of the |//
  //|   msSystem object.                                 |//
  //|                                                    |//
  //|   In order to create a valid kernel the following  |//
  //|   argumens must be provided this method:           |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|		TimeInSeconds = ElapsedTime()                |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

 double	ElapsedTime( void );

 /*/\/\/\/\/\/\/\/\*/
 /*  System Output */
 /*\/\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|				     *  Prompt  *                    |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Outputs to a device a character message contain- |//
  //|   ing delimeters. These delimeters are replaced by |//
  //|   the variable input parameters passed to prompt.  |//
  //|   (Like printf.)                                   |//
  //|                                                    |//
  //|   This method should be altered if a special       |//
  //|   device either than stderr is desired to be used  |//
  //|   as an output prompt.                             |//
  //|                                                    |//
  //|   The arguments to this method are:                |//
  //|                                                    |//
  //|   <* PromptStr *>                                  |//
  //|   A string possibly containing delimeters that     |//
  //|   is to be output to the user.                     |//
  //|                                                    |//
  //|   <* varArgs *>                                    |//
  //|   A variable set of arguments to be placed into    |//
  //|   the prompt string.                               |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|		Prompt(PromptStr, varArgs)                   |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

 void	Prompt(const char*, ...);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|				   *  Progress  *                    |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   This method is called by the mean shift library  |//
  //|   methods during the application of a specific     |//
  //|   algorithm in which the progress of the algorithm |//
  //|   is indicated. Its main use is for a multi-thre-  |//
  //|   aded programming envioronment. Namely, it is     |//
  //|   called fairly frequently by the mean shift       |//
  //|   library methods. It then can be used to grace-   |//
  //|   fully halt the algorithm, or to simply update    |//
  //|   a progress bar.                                  |//
  //|                                                    |//
  //|   This method depends strongly on the interface    |//
  //|   and therefore must be re-implemented to accom-   |//
  //|   odate ones needs.                                |//
  //|                                                    |//
  //|   To facilitate a multi-threaded enviornment       |//
  //|   the prompt function returns a value that         |//
  //|   indicates to the mean shift method whether       |//
  //|   to continue execution. EL_HALT is returned       |//
  //|   when the mean shift procedure is to stop         |//
  //|   execution and EL_OKAY is returned otherwise.     |//
  //|                                                    |//
  //|   The arguments to this method are:                |//
  //|                                                    |//
  //|   <* PercentComplete *>                            |//
  //|   A floating point number that indicates the perc- |//
  //|   ent complete of the algorithm. PercentComplete   |//
  //|   takes a value between zero and one.              |//
  //|                                                    |//
  //|   <* SystemState *>                                |//
  //|   Indicates the system state. It is EL_HALT        |//
  //|   when the mean shift method is to halt execution  |//
  //|   and it is EL_OKAY otherwise.                     |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|	 SystemState = Progress(PercentComplete)         |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

 ErrorLevel Progress(float);

 private:

	 //Timer object...
	 time_t currentTime;

};

#endif
