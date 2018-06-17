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

class Util {

     public:
          static void ParseBranchesByType(RNAStructure *rnaStructBase, RNAStructure::BaseData ** &bdArray, int *btSizes);
          static bool WriteBranchFiles(RNAStructure::BaseData ** &bdArray, int *bdSizes, const RuntimeConfig_t &runtimeConfig);

};

#endif 
