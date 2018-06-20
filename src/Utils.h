/* Utils.h : 
   Utility class with static methods to process data in the program.
   Author:  Maxie D. Schmidt (maxieds@gmail.com)
   Created: 2018.06.17
*/ 

#ifndef __UTILS_H__
#define __UTILS_H__

#include <vector>
using std::vector;

#include "RuntimeConfig.h"
#include "RNAStructure.h"

#define PSPLOT_DIVIDER_STROKE_SIZE     (10)

class Util {

     public:
          static void ParseBranchesByType(RNAStructure *rnaStructBase, 
                                          RNAStructure::BaseData ** &bdArray, int *btSizes);
          static bool WriteBranchFiles(RNAStructure::BaseData ** &bdArray, int *bdSizes, 
                                       const RuntimeConfig_t &runtimeConfig);
          static bool WriteBranchDotBracketFiles(RNAStructure::BaseData ** &bdArray, int *bdSizes, 
                                                 const RuntimeConfig_t &runtimeConfig);
          static bool GenerateDomainPSPlot(const char *outputFile, RuntimeConfig_t runtimeConfig);


     private:
          static void writePSPlotHeaderInfo(FILE *fp);
          static void writePSPlotDomainData(FILE *fp, const char *dotBracketSourceFile); 
          static void writePSPlotActionData(FILE *fp); 
          static void writePSPlotFooterInfo(FILE *fp);  

}; 

#endif 
