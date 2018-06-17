/* BranchTypeIdentification.cpp 
   Author:  Maxie D. Schmidt (maxieds@gmail.com)
   Created: 2018.06.16
*/ 

#include <stdio.h>
#include <stdlib.h>
#include "BranchTypeIdentification.h"


RNABranchType_t::RNABranchType_t(BranchID_t bid = BRANCH_UNDEFINED, class RNAStructure::BaseData *bparent = NULL) {
     branchID = bid;
     branchParent = bparent;
} 
 
RNABranchType_t & RNABranchType_t::operator=(const BranchID_t &rhs) {
     setBranchID(rhs);
     branchParent = NULL;
} 
         
bool RNABranchType_t::operator==(const RNABranchType_t &rhs) const {
     return branchID == rhs.branchID && branchParent == rhs.branchParent;
}
          
bool RNABranchType_t::operator==(const BranchID_t &rhs) const {
     return branchID == rhs;
}

const BranchID_t & RNABranchType_t::getBranchID() const {
     return branchID;
} 
          
void RNABranchType_t::setBranchID(BranchID_t bid) {
     branchID = bid;
}
          
const RNAStructure::BaseData* RNABranchType_t::getBranchParent() const {
     return branchParent;
}
          
void RNABranchType_t::setBranchParent(class RNAStructure::BaseData* bparent) {
     branchParent = bparent;
}

void RNABranchType_t::SetBranchColor(cairo_t * &cr, BranchID_t bt) {
     switch(bt) {
          case BRANCH_UNDEFINED:
               cairo_set_source_rgb(cr, 0.0, 0.0, 0.0); 
               break;
          case BRANCH1:
               cairo_set_source_rgb(cr, 92.0 / 255, 160.0 / 255, 215.0 / 255);
               break;
          case BRANCH2:
               cairo_set_source_rgb(cr, 183.0 / 255, 127.0 / 255, 77.0 / 255);
               break;
          case BRANCH3:
               cairo_set_source_rgb(cr, 243.0 / 255, 153.0 / 255, 193.0 / 255);
               break;
          case BRANCH4:
               cairo_set_source_rgb(cr, 123.0 / 255, 204.0 / 255, 153.0 / 255);
               break;
          default:
               break;
     }
}

bool RNABranchType_t::PerformBranchClassification(class RNAStructure *rnaStructBase, unsigned int alength) {

     if(alength < 4)
          return false;     

     // first we determine the four most enclosing arcs on the circle: 
     RNAStructure::BaseData* mostEnclosingArcs[4] = {NULL, NULL, NULL, NULL};
     unsigned int mostEnclosingArcsSize = 0;
     for(int rs = 0; rs < alength; rs++) {
          //fprintf(stderr, "[1] rs=%d\n", rs); 
          RNAStructure::BaseData* rnaStruct = rnaStructBase->GetBaseAt(rs);
          if(rnaStruct->m_pair == RNAStructure::UNPAIRED) 
               continue;
          else if(mostEnclosingArcsSize == 0) {
               mostEnclosingArcs[0] = rnaStructBase->GetBaseAt(rs);
               //fprintf(stderr, "mostEnclosingArcs[0]=%p\n", mostEnclosingArcs[0]);
               mostEnclosingArcsSize++;
               continue;
          }
          bool isEnclosedInLargerBranch = false;
          for(int mea = 0; mea < mostEnclosingArcsSize; mea++) {
               if(rnaStruct->isContainedIn(*(mostEnclosingArcs[mea]))) {
                    isEnclosedInLargerBranch = true;
                    rnaStructBase->GetBranchTypeAt(rs).setBranchID((BranchID_t) (mea + 1));
                    rnaStructBase->GetBranchTypeAt(rs).setBranchParent(mostEnclosingArcs[mea]);
                    break;
               }
          }
          if(isEnclosedInLargerBranch) 
               continue; // this cannot be an outer containing branch
          RNAStructure::BaseData *currentBaseData = rnaStructBase->GetBaseAt(rs);;
          unsigned int pairDistance = ABS(MAX(currentBaseData->m_pair, currentBaseData->m_index) - 
                                          MIN(currentBaseData->m_pair, currentBaseData->m_index));
          for(int mea = 0; mea < mostEnclosingArcsSize; mea++) { 
               RNAStructure::BaseData *meaBaseData = mostEnclosingArcs[mea];
               unsigned int meaPairDistance = ABS(MAX(meaBaseData->m_pair, meaBaseData->m_index) - 
                                                  MIN(meaBaseData->m_pair, meaBaseData->m_index)); 
               bool needToResort = false;
               if(meaPairDistance < pairDistance && meaBaseData->m_pair > meaBaseData->m_index) {
                    //fprintf(stderr, "[2] meaPD=%d, PD=%d\n", meaPairDistance, pairDistance);
                    if(mostEnclosingArcsSize < 4) {
                         mostEnclosingArcs[mostEnclosingArcsSize] = mostEnclosingArcs[mea];
                         mostEnclosingArcsSize++;
                         needToResort = true;
                    }
                    mostEnclosingArcs[mea] = rnaStructBase->GetBaseAt(rs);
                    rnaStructBase->GetBranchTypeAt(rs).branchID = (BranchID_t) (mea + 1);
                    if(needToResort) {
                         qsort(&mostEnclosingArcs[0], mostEnclosingArcsSize, sizeof(RNAStructure::BaseData*), compareMostEnclosingArcs);
                    }
                    //for(int i = 0; i < mostEnclosingArcsSize; i++) {
                    //     fprintf(stderr, "   => %d, %d, %p\n", i, needToResort, mostEnclosingArcs[i]);
                    //     fprintf(stderr, "   => %d: m_idx=%d, mPD=%d\n", i, mostEnclosingArcs[i]->m_index, 
                    //             mostEnclosingArcs[i]->getPairDistance());
                    //}
                    break;
               }
          }
     }
     // now that we've identified most of the the enclosing branches, 
     // we reset the branch types by number on all (except for the nubbins, 
     // see below) entries in the array: 
     for(int rs = 0; rs < alength; rs++) {
          //if(mostEnclosingArcsSize < 4) {
          //     rnaStructBase->GetBranchTypeAt(rs).setBranchID(BRANCH_UNDEFINED);
          //     rnaStructBase->GetBranchTypeAt(rs).setBranchParent(NULL);
          //     continue;
          //}
          bool isNubbin = true;
          for(int mea = 0; mea < mostEnclosingArcsSize; mea++) { 
               if(rnaStructBase->GetBaseAt(rs)->isContainedIn(*(mostEnclosingArcs[mea]))) {
                    rnaStructBase->GetBranchTypeAt(rs).setBranchID((BranchID_t) (mea + 1));
                    rnaStructBase->GetBranchTypeAt(rs).setBranchParent(mostEnclosingArcs[mea]);
                    isNubbin = false;
                    break;
               }
          }
          if(isNubbin) {
               unsigned int localIndex = rnaStructBase->GetBaseAt(rs)->m_index;
               unsigned int localPair = rnaStructBase->GetBaseAt(rs)->m_pair;
               if(localPair == RNAStructure::UNPAIRED)
                    localPair = localIndex;
               IntIndexPair_t closestEnclosingArcs[4];
               for(int mea = mostEnclosingArcsSize - 1; mea >= 0; mea--) {
                    IntIndexPair_t iip;
                    iip.index = mea;
                    int dist1 = ABS(localIndex - mostEnclosingArcs[mea]->m_pair);
                    int dist2 = ABS(mostEnclosingArcs[mea]->m_index - localPair);
                    int dist3 = ABS(mostEnclosingArcs[mea]->m_pair - localIndex);
                    iip.intValue = MIN3(dist1, dist2, dist3);
                    closestEnclosingArcs[mea] = iip;
               }
               qsort(closestEnclosingArcs, mostEnclosingArcsSize, sizeof(IntIndexPair_t), compareIntegerIndexPair);
               rnaStructBase->GetBranchTypeAt(rs).setBranchID((BranchID_t) (closestEnclosingArcs[3].index + 1));
               rnaStructBase->GetBranchTypeAt(rs).setBranchParent(mostEnclosingArcs[closestEnclosingArcs[3].index]);
          }
     }
     //for(int i = 0; i < mostEnclosingArcsSize; i++) {
     //     fprintf(stderr, "   => %d: m_idx=%u, mPD=%u\n", i, mostEnclosingArcs[i]->m_index, 
     //             mostEnclosingArcs[i]->getPairDistance());
     //}
     if(mostEnclosingArcsSize < 4)
          return false;
     return true;

} 
