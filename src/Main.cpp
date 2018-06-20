/* Main.cpp : 
   Main runner for the RNABranchIDViz program. 
   Author:  Maxie D. Schmidt
   Created: 2018.06.17
*/ 

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#include "RuntimeConfig.h"
#include "BranchTypeIdentification.h"
#include "DiagramImage.h"
#include "RNAStructure.h"
#include "Utils.h"

int main(int argc, char **argv) {

     RuntimeConfig_t runtimeOptions;
     if(!runtimeOptions.parseRuntimeArgs(argc, argv)) { 
          RuntimeConfig_t::PrintUsage(argv[0]);
          exit(-1);
     }

     RNAStructure *rnaStructBase = RNAStructure::CreateFromFile(runtimeOptions.getBaseFilePathOption(), false); 
     RNAStructure::BaseData **bdArray = (RNAStructure::BaseData**) malloc(sizeof(RNAStructure::BaseData*) * 4); 
     int *bdSizes = (int *) malloc(sizeof(int) * 4);
     for(int bd = 0; bd < NUM_BRANCHES; bd++) { 
          bdSizes[bd] = 0;
          bdArray[bd] = (RNAStructure::BaseData*) malloc(sizeof(RNAStructure::BaseData) * rnaStructBase->GetLength()); 
     }
     Util::ParseBranchesByType(rnaStructBase, bdArray, bdSizes);
     Util::WriteBranchFiles(bdArray, bdSizes, runtimeOptions); 
     Util::WriteBranchDotBracketFiles(bdArray, bdSizes, runtimeOptions); 

     if(runtimeOptions.getOutputImagesOption()) { 
          DiagramImage_t *rnaSVDiagramImage = new DiagramImage_t(rnaStructBase);
          char rsvImageOutFile[MAX_FILEPATH_LENGTH]; 
          snprintf(rsvImageOutFile, MAX_FILEPATH_LENGTH, "%s-RNAStructViz.png", runtimeOptions.getBaseFilePathNoCTOption());
          rnaSVDiagramImage->writePNGImage(rsvImageOutFile); 
          delete rnaSVDiagramImage;
          char viennaPostscriptOutFile[MAX_FILEPATH_LENGTH];
          snprintf(viennaPostscriptOutFile, MAX_FILEPATH_LENGTH, "%s-ViennaSubdomainPlot.ps",runtimeOptions.getBaseFilePathNoCTOption());
          Util::GenerateDomainPSPlot(viennaPostscriptOutFile, runtimeOptions);
     } 

     free(rnaStructBase);
     free(bdSizes); 
     for(int bd = 0; bd < NUM_BRANCHES; bd++) { 
          free(bdArray[bd]);
     } 
     free(bdArray);

     return 0;

}
