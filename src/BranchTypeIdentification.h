/* BranchTypeIdentification.h : 
   A helper class to identify which of the four branch types 
   to which the RNAStructure and/or pairing belongs.
   Author:  Maxie D. Schmidt (maxieds@gmail.com)
   Created: 2018.06.16
*/

#ifndef __RNABRANCHTYPEIDENT_H__
#define __RNABRANCHTYPEIDENT_H__

#include <FL/fl_draw.H>
#include "RNAStructure.h" 

#define FL_COLOR_BRANCH_UNDEFINED     (FL_BLACK)
#define FL_COLOR_BRANCH1              (fl_rgb_color(92, 160, 215))
#define FL_COLOR_BRANCH2              (fl_rgb_color(183, 127, 77))
#define FL_COLOR_BRANCH3              (fl_rgb_color(243, 153, 193))
#define FL_COLOR_BRANCH4              (fl_rgb_color(123, 204, 153))

#define NUM_BRANCHES                  (4)

typedef enum {
     BRANCH1 = 1, 
     BRANCH2 = 2, 
     BRANCH3 = 3, 
     BRANCH4 = 4, 
     BRANCH_UNDEFINED = 0
} BranchID_t;

typedef struct {
     int index;
     int intValue;
} IntIndexPair_t;

inline int compareIntegerIndexPair(const void *ip1, const void *ip2) {
     int i1 = ((IntIndexPair_t*) ip1)->intValue, i2 = ((IntIndexPair_t*) ip2)->intValue;
     if(i1 < i2) 
          return -1;
     else if(i1 == i2)
          return 0;
     else
          return 1;
}

inline int compareMostEnclosingArcs(const void *arc1, const void *arc2) {
     RNAStructure::BaseData *rnaStruct1 = (RNAStructure::BaseData *) arc1;
     RNAStructure::BaseData *rnaStruct2 = (RNAStructure::BaseData *) arc2;
     unsigned int arc1PairDist = rnaStruct1->getPairDistance();
     unsigned int arc2PairDist = rnaStruct2->getPairDistance();
     if(arc1PairDist < arc2PairDist)
          return -1;
     else if(arc1PairDist == arc2PairDist)
          return 0;
     else
          return 1;
}

class RNABranchType_t {

     public: 
          RNABranchType_t(BranchID_t bid, class RNAStructure::BaseData *bparent); 
 
          RNABranchType_t & operator=(const BranchID_t &rhs); 
         
          bool operator==(const RNABranchType_t &rhs) const;
          bool operator==(const BranchID_t &rhs) const;

          const BranchID_t & getBranchID() const;
          void setBranchID(BranchID_t bid); 
          
          const RNAStructure::BaseData* getBranchParent() const;
          void setBranchParent(class RNAStructure::BaseData* bparent); 

          static void SetFLBranchColor(BranchID_t bt);

          static bool PerformBranchClassification(class RNAStructure *rnaStructArray, unsigned int alength);

     protected:
          BranchID_t branchID;
          RNAStructure::BaseData *branchParent;

};

#endif
