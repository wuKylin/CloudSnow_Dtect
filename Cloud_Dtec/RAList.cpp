// RAList.cpp: implementation of the RAList class.
//
//////////////////////////////////////////////////////////////////////


#include "RAList.h"
//include needed libraries
#include	<stdio.h>
#include	<assert.h>
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
/*Constructs a RAList object.                          */
/*******************************************************/
/*Post:                                                */
/*      - a RAlist object has been properly constru-   */
/*        cted.                                        */
/*******************************************************/

RAList::RAList( void )
{
	//initialize label and link
	label			= -1;
	next			= NULL;

	//initialize edge strenght weight and count
	edgeStrength	= 0; // 边缘强处
	edgePixelCount	= 0; // 边缘像素点的个数
}

/*******************************************************/
/*Class Destructor                                     */
/*******************************************************/
/*Destructrs a RAList object.                          */
/*******************************************************/
/*Post:                                                */
/*      - the RAList object has been properly dest-    */
/*        ructed.                                      */
/*******************************************************/

RAList::~RAList( void )
{
	//do nothing
}

/*******************************************************/
/*Insert                                               */
/*******************************************************/
/*Insert a region node into the region adjacency list. */
/*******************************************************/
/*Pre:                                                 */
/*      - entry is a node representing a connected re- */
/*        gion                                         */
/*Post:                                                */
/*      - entry has been inserted into the region adj- */
/*        acency list if it does not already exist     */
/*        there.                                       */
/*      - if the entry already exists in the region    */
/*        adjacency list 1 is returned otherwise 0 is  */
/*        returned.                                    */
/*******************************************************/

// 其就是根据label的一个有序链表，label从小到大
int RAList::Insert(RAList *entry)
{

	//if the list contains only one element
	//then insert this element into next
	if(!next)
	{
		//insert entry
		next		= entry;  // 
		entry->next = NULL;

		//done
		return 0;
	}

	//traverse the list until either:

	//(a) entry's label already exists - do nothing
	//(b) the list ends or the current label is
	//    greater than entry's label, thus insert the entry
	//    at this location


	//check first entry, 头插法，RAList的next是指向头结点
	if(next->label > entry->label)
	{
		//insert entry into the list at this location
		entry->next	= next;
		next		= entry; 

		//done
		return 0;
	}

	// next 存储一个头结点
	//check the rest of the list...
	exists	= 0;
	cur		= next; 
	while(cur)
	{
		if(entry->label == cur->label)
		{
			//node already exists
			exists = 1; // 该元素已经存在
			break;
		}
		else if((!(cur->next))||(cur->next->label > entry->label))
		{
			//insert entry into the list at this location
			entry->next	= cur->next;//
			cur->next	= entry;//插入到合适的位置
			break;
		}

		//traverse the region adjacency list
		cur = cur->next;// 
	}

	//done. Return exists indicating whether or not a new node was
	//      actually inserted into the region adjacency list.
	return (int)(exists);

}
