/* RuntimeConfig.cpp 
   Author:  Maxie D. Schmidt (maxieds@gmail.com)
   Created: 2018.06.17
*/ 

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>

#include "RuntimeConfig.h"
#include "BranchTypeIdentification.h"

RuntimeConfig_t::RuntimeConfig_t() : quiet(0), debug(0), outputImages(1), 
     outputFasta(0), renumberCTIndices(1) {
     baseFilePath[0] = '\0';
     baseFilePathNoCT[0] = '\0';
} 

int RuntimeConfig_t::getQuietOption() const {
     return quiet;
} 

int RuntimeConfig_t::getDebugOption() const {
     return debug;
} 

int RuntimeConfig_t::getOutputFASTAOption() const {
     return outputFasta;
}

int RuntimeConfig_t::getRenumberOption() const {
     return renumberCTIndices;
} 

int RuntimeConfig_t::getOutputImagesOption() const {
     return outputImages;
}

const char * RuntimeConfig_t::getBaseFilePathOption() const {
     return baseFilePath;
}

const char * RuntimeConfig_t::getBaseFilePathNoCTOption() const {
     return baseFilePathNoCT;
}
bool RuntimeConfig_t::parseRuntimeArgs(int argc, char **argv) {
     
     if(argc < 2 || strlen(argv[1]) < 3)
          return false;
     strncpy(baseFilePath, argv[1], MAX_FILEPATH_LENGTH);      
     strncpy(baseFilePathNoCT, argv[1], MAX_FILEPATH_LENGTH);
     baseFilePathNoCT[strlen(argv[1]) - 3] = '\0';
     argc = argc - 1;
     argv[1] = argv[0];
     argv = &argv[1];

     static struct option longOptions[] = {
          {"quiet", no_argument, (int *) &quiet, 1}, 
          {"debug", no_argument, (int *) &debug, 1}, 
          {"output-fasta", no_argument, (int *) &outputFasta, 1}, 
          {"no-renumber-CT", no_argument, (int *) &renumberCTIndices, 0}, 
          {"no-images", no_argument, (int *) &outputImages, 0}, 
          {0, 0, 0, 0}
     }; 
     int optc, optionIndex = 0;
     while(true) {
          optc = getopt_long(argc, argv, ":", longOptions, &optionIndex); 
          if(optc == -1)
               break;
          switch(optc) { 
               case 0: // option sets a flag, do nothing else now
                    break;
               default:
                    return false;
          }
     } 
     if(optind < argc)
          return false; // unknown commandline options
     return true;
}

char * RuntimeConfig_t::getBranchTypeOutputFile(BranchID_t bid, const char *fileExt = "ct") const { 
     if(bid == BRANCH_UNDEFINED)
          return NULL;
     char *outputFilePath = (char *) malloc(sizeof(char) * MAX_FILEPATH_LENGTH); 
     snprintf(outputFilePath, MAX_FILEPATH_LENGTH, "%s-branch%02d.%s", baseFilePathNoCT, bid, fileExt);
     return outputFilePath; 
}

void RuntimeConfig_t::PrintUsage(char *progName) {
     fprintf(stderr, "Usage: %s CTFileName [--quiet] [--debug] [--output-fasta] [--no-renumber-CT] [--no-images]\n", progName);
} 


