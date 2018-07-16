/* Utils.h : 
   Utility class with static methods to process data in the program.
   Author:  Maxie D. Schmidt (maxieds@gmail.com)
   Created: 2018.06.17
*/ 

#ifndef __UTILS_H__
#define __UTILS_H__

#include <vector>
#include <algorithm>
using std::vector;
using std::sort;

#include "RuntimeConfig.h"
#include "RNAStructure.h"

#define PSPLOT_DIVIDER_STROKE_SIZE     (10)

class Util {

     public:
          static float getDomainRed(BranchID_t bid);
          static float getDomainGreen(BranchID_t bid);
          static float getDomainBlue(BranchID_t bid);

          static vector<RNAStructure::BaseData *> getEnclosingArcs(RNAStructure * &rnaStructBase, bool removeTopFour);
          static void ParseBranchesByType(RNAStructure * &rnaStructBase, 
                                          RNAStructure::BaseData ** &bdArray, int *btSizes, RuntimeConfig_t runtimeConfig);
          static bool WriteBranchFiles(RNAStructure::BaseData ** &bdArray, int *bdSizes, 
                                       const RuntimeConfig_t &runtimeConfig);
          static bool WriteFullDotBracketFile(RNAStructure *rnaStructBase, RuntimeConfig_t runtimeConfig); 
          static bool WriteBranchDotBracketFiles(RNAStructure::BaseData ** &bdArray, int *bdSizes, 
                                                 const RuntimeConfig_t &runtimeConfig);
          static bool WriteBranchFASTAFiles(RNAStructure::BaseData ** &bdArray, int *bdSizes, 
                                            const RuntimeConfig_t &runtimeConfig); 
          static bool GenerateDomainPSPlot(const char *outputFile, RuntimeConfig_t runtimeConfig);


     private:
          static void writePSPlotHeaderInfo(FILE *fp);
          static void writePSPlotDomainData(FILE *fp, const char *dotBracketSourceFile); 
          static void writePSPlotActionData(FILE *fp); 
          static void writePSPlotFooterInfo(FILE *fp);  

     public:
          typedef struct {
               inline bool operator()(RNAStructure::BaseData *bd1, RNAStructure::BaseData *bd2) { // sorts in decreasing order (i.e., largest arcs first): 
                    int bd1ArcDist = (bd1->m_pair == RNAStructure::UNPAIRED) ? 0 : MAX(bd1->m_index, bd1->m_pair) - MIN(bd1->m_index, bd1->m_pair);
                    int bd2ArcDist = (bd2->m_pair == RNAStructure::UNPAIRED) ? 0 : MAX(bd2->m_index, bd2->m_pair) - MIN(bd2->m_index, bd2->m_pair);
                    return bd1ArcDist > bd2ArcDist;
               }
          } RNAStructureBaseDataArcLengthSort;

          typedef struct {
               inline bool operator()(RNAStructure::BaseData *bd1, RNAStructure::BaseData *bd2) { 
                    return bd1->m_index < bd2->m_index;
               }
          } RNAStructureBaseDataIndexSort;

}; 

#endif 
