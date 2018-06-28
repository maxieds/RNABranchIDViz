/* RuntimeConfig.h : 
   Defines runtime parameters for the RNABranchIDViz program. 
   Author:  Maxie D. Schmidt (maxieds@gmail.com)
   Created: 2018.06.17
*/

#ifndef __RUNTIME_CONFIG_H__
#define __RUNTIME_CONFIG_H__

#include "BranchTypeIdentification.h"

#define MAX_FILEPATH_LENGTH      (256)

typedef enum {
     FLAT, 
     FLAT_UNCIRCLED, 
     CIRCLE, 
     CIRCLE_UNCIRCLED,
     DEFAULT
} DrawImageType_t; 

class RuntimeConfig_t {

     private:
          int quiet;
          int debug;
          int renumberCTIndices;
          int outputImages;
          char baseFilePath[MAX_FILEPATH_LENGTH];
          char baseFilePathNoCT[MAX_FILEPATH_LENGTH];

     public:
          RuntimeConfig_t();
          
          bool getQuietOption() const;
          bool getDebugOption() const;
          bool getRenumberOption() const;
          bool getOutputImagesOption() const;
          const char * getBaseFilePathOption() const;
          const char * getBaseFilePathNoCTOption() const;

          bool parseRuntimeArgs(int argc, char **argv); 
          char * getBranchTypeOutputFile(BranchID_t bid, const char *fileExt) const;

          static void PrintUsage(char *progName);

};


#endif
