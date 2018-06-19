/* Utils.cpp 
   Author:  Maxie D. Schmidt (maxieds@gmail.com)
   Created: 2018.06.17
*/ 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#include <map>
#include <utility>
using std::map;
using std::pair;

#include "Utils.h"

void Util::ParseBranchesByType(RNAStructure *rnaStructBase, RNAStructure::BaseData ** &bdArray, int *btSizes) {
     for(int rs = 0; rs < rnaStructBase->GetLength(); rs++) { 
          RNAStructure::BaseData curBaseData = *(rnaStructBase->GetBaseAt(rs));
          fprintf(stderr, "Base @ %d: %u -> %u\n", rs, curBaseData.m_index, curBaseData.m_pair); 
          switch(rnaStructBase->GetBranchTypeAt(rs).getBranchID()) { 
               case BRANCH1:
                    bdArray[0][btSizes[0]] = curBaseData;
                    btSizes[0]++;
                    break;
               case BRANCH2:
                    bdArray[1][btSizes[1]] = curBaseData;
                    btSizes[1]++;
                    break;
               case BRANCH3:
                    bdArray[2][btSizes[2]] = curBaseData;
                    btSizes[2]++;
                    break;
               case BRANCH4:
                    bdArray[3][btSizes[3]] = curBaseData;
                    btSizes[3]++;
                    break;
               default:
                    break;
          }
     }
}

bool Util::WriteBranchFiles(RNAStructure::BaseData ** &bdArray, int *bdSizes, const RuntimeConfig_t &runtimeConfig) {
     for(int bt = 0; bt < NUM_BRANCHES; bt++) { 
          // generate pair mappings (i.e., renumberings within the local file): 
          map <unsigned int, unsigned int> bdPairingsMap;
          bdPairingsMap.insert(pair <unsigned int, unsigned int>(RNAStructure::UNPAIRED, 0));
          for(int bd = 0; bd < bdSizes[bt]; bd++) { 
               //fprintf(stderr, "Mapping %u -> %d [paired with %u]\n", bdArray[bt][bd].m_index, bd + 1, bdArray[bt][bd].m_pair);
               bdPairingsMap.insert(pair <unsigned int, unsigned int>(bdArray[bt][bd].m_index, bd + 1)); 
          } 
          char *branchOutputFile = runtimeConfig.getBranchTypeOutputFile((BranchID_t) (bt + 1));
          if(!runtimeConfig.getQuietOption()) { 
               fprintf(stdout, "BRANCH TYPE #%d INCLUSIONS [\"%s\"]:\n", bt + 1, branchOutputFile); 
          } 
          FILE *btfp = fopen(branchOutputFile, "w+"); 
          if(!btfp) {
               perror("fopen"); 
               free(branchOutputFile);
               return false;
          }
          // write file header information: 
          char writeBuf[MAX_FILEPATH_LENGTH];
          snprintf(writeBuf, MAX_FILEPATH_LENGTH, "Filename: %s\n", branchOutputFile); 
          fwrite(writeBuf, sizeof(char), strlen(writeBuf), btfp); 
          snprintf(writeBuf, MAX_FILEPATH_LENGTH, "Branch #%d CT file written by RNABranchIDViz\n", bt + 1);
          fwrite(writeBuf, sizeof(char), strlen(writeBuf), btfp); 
          snprintf(writeBuf, MAX_FILEPATH_LENGTH, " %d   dG =     0.00  [initially     0.0]\n", bdSizes[bt]);
          fwrite(writeBuf, sizeof(char), strlen(writeBuf), btfp); 
          // write pair information: 
          for(int bd = 0; bd < bdSizes[bt]; bd++) { 
               int curIndex = runtimeConfig.getRenumberOption() ? bd + 1 : bdArray[bt][bd].m_index; 
               int pairIndex = runtimeConfig.getRenumberOption() ? bdPairingsMap.at(bdArray[bt][bd].m_pair) : \
                               bdArray[bt][bd].m_pair; 
               const char *newline = (bd + 1 == bdSizes[bt]) ? "" : "\n";
               snprintf(writeBuf, MAX_FILEPATH_LENGTH, "% 5d %c % 5d % 5d % 5d % 5d%s", 
                        curIndex, bdArray[bt][bd].getBaseChar(), curIndex - 1, curIndex + 1, pairIndex, curIndex, newline); 
               fwrite(writeBuf, sizeof(char), strlen(writeBuf), btfp); 
               if(!runtimeConfig.getQuietOption()) { 
                    curIndex = bdArray[bt][bd].m_index + 1; 
                    snprintf(writeBuf, MAX_FILEPATH_LENGTH, "% 5d %c % 5d % 5d % 5d % 5d%s", 
                        curIndex, bdArray[bt][bd].getBaseChar(), curIndex - 1, curIndex + 1, bdArray[bt][bd].m_pair + 1, curIndex, newline); 
                    fprintf(stdout, "%s", writeBuf); 
               }
          } 
          if(!runtimeConfig.getQuietOption()) { 
               fprintf(stdout, "\n"); 
          }
          free(branchOutputFile); 
     } 
     return true;
}

bool Util::GenerateBranchDrawImages(const char *outPrefix, DrawImageType_t drawSpec) { 
     return true;
}

