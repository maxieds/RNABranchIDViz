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

float Util::getDomainRed(BranchID_t bid) { 
     switch(bid) { 
          case BRANCH1:
               return 92.0 / 255;
          case BRANCH2: 
               return 183.0 / 255;
          case BRANCH3:
               return 243.0 / 255;
          case BRANCH4:
               return 123.0 / 255;
          default:
               break;
     }
     return 0.0;
}

float Util::getDomainGreen(BranchID_t bid) { 
     switch(bid) { 
          case BRANCH1:
               return 160.0 / 255;
          case BRANCH2: 
               return 127.0 / 255;
          case BRANCH3:
               return 153.0 / 255;
          case BRANCH4:
               return 204.0 / 255;
          default:
               break;
     }
     return 0.0;
}

float Util::getDomainBlue(BranchID_t bid) { 
     switch(bid) { 
          case BRANCH1:
               return 215.0 / 255;
          case BRANCH2: 
               return 77.0 / 255;
          case BRANCH3:
               return 193.0 / 255;
          case BRANCH4:
               return 153.0 / 255;
          default:
               break;
     }
     return 0.0;
}

vector<RNAStructure::BaseData *> Util::getEnclosingArcs(RNAStructure * &rnaStructBase, bool removeTopFour = false) { 
     vector<RNAStructure::BaseData *> rBaseDataVec;
     for(int bd = 0; bd < rnaStructBase->GetLength(); bd++) {
          bool isAlreadyEnclosed = false;
          RNAStructure::BaseData *curBaseData = rnaStructBase->GetBaseAt(bd);
          for(int v = 0; v < rBaseDataVec.size(); v++) {
               if(curBaseData->m_pair == RNAStructure::UNPAIRED) { 
                    isAlreadyEnclosed = true;
                    break;
               }
               else if(curBaseData->isContainedIn(*(rBaseDataVec[v]))) {
                    isAlreadyEnclosed = true;
                    break;
               }
               else if(rBaseDataVec[v]->isContainedIn(*curBaseData)) {
                    rBaseDataVec[v] = curBaseData;
                    isAlreadyEnclosed = true;
                    break;
               }
          }
          if(!isAlreadyEnclosed && (rBaseDataVec.size() > 0 || curBaseData->m_pair != RNAStructure::UNPAIRED)) {
               rBaseDataVec.push_back(curBaseData);
          }
     }
     sort(rBaseDataVec.begin(), rBaseDataVec.end(), RNAStructureBaseDataArcLengthSort());
     if(removeTopFour) {
          rBaseDataVec.erase(rBaseDataVec.begin(), rBaseDataVec.begin() + 4);
     }
     return rBaseDataVec;
}


void Util::ParseBranchesByType(RNAStructure * &rnaStructBase, RNAStructure::BaseData ** &bdArray, int *btSizes, RuntimeConfig_t runtimeConfig) {
     for(int rs = 0; rs < rnaStructBase->GetLength(); rs++) { 
          RNAStructure::BaseData curBaseData = *(rnaStructBase->GetBaseAt(rs));
          if(runtimeConfig.getDebugOption()) { 
               fprintf(stderr, "Base @ %d: %u -> %u\n", rs, curBaseData.m_index, curBaseData.m_pair); 
          }
          switch(rnaStructBase->GetBranchTypeAt(rs)->getBranchID()) { 
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
               bdPairingsMap.insert(pair <unsigned int, unsigned int>(bdArray[bt][bd].m_index, bd + 1)); 
          } 
          char *branchOutputFile = runtimeConfig.getBranchTypeOutputFile((BranchID_t) (bt + 1), "ct");
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
          fclose(btfp);
          free(branchOutputFile); 
     } 
     return true;
}

bool Util::WriteFullDotBracketFile(RNAStructure *rnaStructBase, RuntimeConfig_t runtimeConfig) {
     char dbOutputFile[MAX_FILEPATH_LENGTH];
     snprintf(dbOutputFile, MAX_FILEPATH_LENGTH, "%s.dot", runtimeConfig.getBaseFilePathNoCTOption());
     FILE *outfp = fopen(dbOutputFile, "w+"); 
     if(!outfp) { 
          fprintf(stderr, "Unable to open \"%s\" for writing : %s\n", dbOutputFile, strerror(errno));
          return false;
     }
     char writeBuf[MAX_FILEPATH_LENGTH];
     snprintf(writeBuf, MAX_FILEPATH_LENGTH, "> Full DOT file for %s\n", runtimeConfig.getBaseFilePathOption());
     fwrite(writeBuf, sizeof(char), strlen(writeBuf), outfp); 
     for(int b = 0; b < rnaStructBase->GetLength(); b++) { 
          const char *newline = (b + 1) == rnaStructBase->GetLength() ? "\n" : "";
          snprintf(writeBuf, MAX_FILEPATH_LENGTH, "%c%s", rnaStructBase->GetBaseAt(b)->getBaseChar(), newline);
          fwrite(writeBuf, sizeof(char), strlen(writeBuf), outfp);
     }
     for(int b = 0; b < rnaStructBase->GetLength(); b++) { 
          const char *newline = (b + 1) == rnaStructBase->GetLength() ? "\n" : "";
          char pairChar = '.';
          if(rnaStructBase->GetBaseAt(b)->m_pair != RNAStructure::UNPAIRED) { 
               pairChar = (rnaStructBase->GetBaseAt(b)->m_index < rnaStructBase->GetBaseAt(b)->m_pair) ? '(' : ')';
          }
          snprintf(writeBuf, MAX_FILEPATH_LENGTH, "%c%s", pairChar, newline);
          fwrite(writeBuf, sizeof(char), strlen(writeBuf), outfp);
     }
     return true;
}

bool Util::WriteBranchDotBracketFiles(RNAStructure::BaseData ** &bdArray, int *bdSizes, 
                                      const RuntimeConfig_t &runtimeConfig) {
     for(int bt = 0; bt < NUM_BRANCHES; bt++) { 
          // write file header information: 
          char *branchOutputFile = runtimeConfig.getBranchTypeOutputFile((BranchID_t) (bt + 1), "dot");
          FILE *btfp = fopen(branchOutputFile, "w+"); 
          if(!btfp) {
               perror("fopen"); 
               free(branchOutputFile);
               return false;
          }
          char writeBuf[MAX_FILEPATH_LENGTH];
          snprintf(writeBuf, MAX_FILEPATH_LENGTH, ">Filename: %s\n", branchOutputFile); 
          fwrite(writeBuf, sizeof(char), strlen(writeBuf), btfp); 
          // write pair information: 
          for(int bd = 0; bd < bdSizes[bt]; bd++) { 
               const char *newline = (bd + 1 == bdSizes[bt]) ? "\n" : "";
               snprintf(writeBuf, MAX_FILEPATH_LENGTH, "%c%s", bdArray[bt][bd].getBaseChar(), newline);
               fwrite(writeBuf, sizeof(char), strlen(writeBuf), btfp); 
          }
          for(int bd = 0; bd < bdSizes[bt]; bd++) { 
               const char *newline = (bd + 1 == bdSizes[bt]) ? "\n" : "";
               char pairingChar = '.';
               if(bdArray[bt][bd].m_pair != RNAStructure::UNPAIRED) { 
                    pairingChar = (bdArray[bt][bd].m_index < bdArray[bt][bd].m_pair) ? '(' : ')';
               }
               snprintf(writeBuf, MAX_FILEPATH_LENGTH, "%c%s", pairingChar, newline);
               fwrite(writeBuf, sizeof(char), strlen(writeBuf), btfp); 
          }
          fclose(btfp);
          free(branchOutputFile); 
     } 
     return true;
}

bool Util::WriteBranchFASTAFiles(RNAStructure::BaseData ** &bdArray, int *bdSizes, 
                                 const RuntimeConfig_t &runtimeConfig) {
     for(int bt = 0; bt < NUM_BRANCHES; bt++) { 
          // write file header information: 
          char *branchOutputFile = runtimeConfig.getBranchTypeOutputFile((BranchID_t) (bt + 1), "fasta");
          FILE *btfp = fopen(branchOutputFile, "w+"); 
          if(!btfp) {
               perror("fopen"); 
               free(branchOutputFile);
               return false;
          }
          char writeBuf[MAX_FILEPATH_LENGTH];
          snprintf(writeBuf, MAX_FILEPATH_LENGTH, ">Filename: %s\n", branchOutputFile); 
          fwrite(writeBuf, sizeof(char), strlen(writeBuf), btfp); 
          // write pair information: 
          for(int bd = 0; bd < bdSizes[bt]; bd++) { 
               const char *newline = (bd + 1 == bdSizes[bt]) ? "\n" : "";
               snprintf(writeBuf, MAX_FILEPATH_LENGTH, "%c%s", bdArray[bt][bd].getBaseChar(), newline);
               fwrite(writeBuf, sizeof(char), strlen(writeBuf), btfp); 
          }
          fclose(btfp);
          free(branchOutputFile); 
     } 
     return true;
}

