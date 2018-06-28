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

extern "C" {
     #include <ViennaRNA/model.h>
     #include <ViennaRNA/plotting/layouts.h>
     #include <ViennaRNA/utils/basic.h>
     #include <ViennaRNA/utils/alignments.h>
     #include <ViennaRNA/model.h>
     #include <ViennaRNA/fold_vars.h>
     #include <ViennaRNA/gquad.h>
     #include <ViennaRNA/plotting/layouts.h>
     #include <ViennaRNA/alphabet.h>
     #include <ViennaRNA/plotting/structures.h>
     #include <ViennaRNA/io/file_formats.h>
} 

#include "Utils.h"

extern RNAStructure *rnaStructBaseOrig = NULL;

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
     rnaStructBaseOrig = rnaStructBase;
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
               if(bt == 0 && bd == 0) {
               fprintf(stderr, "bt=%d; bd=%d; m_pair=%d; m_index=%d\n", bt, bd, bdArray[bt][bd].m_pair, bdArray[bt][bd].m_index);
               fprintf(stderr, "bt=%d; bd=%d; m_pair=%d; m_index=%d; BTYPE=%d\n", bt, bd, rnaStructBaseOrig->GetBaseAt(783)->m_pair, rnaStructBaseOrig->GetBaseAt(783)->m_index, rnaStructBaseOrig->GetBranchTypeAt(783)->getBranchID());
               fprintf(stderr, "bt=%d; bd=%d; m_pair=%d; m_index=%d; BTYPE=%d\n", bt, bd, rnaStructBaseOrig->GetBaseAt(756)->m_pair, rnaStructBaseOrig->GetBaseAt(756)->m_index, rnaStructBaseOrig->GetBranchTypeAt(756)->getBranchID());}

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
          snprintf(writeBuf, MAX_FILEPATH_LENGTH, "> Filename: %s\n", branchOutputFile); 
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

bool Util::GenerateDomainPSPlot(const char *outputFile, RuntimeConfig_t runtimeConfig) {
     FILE *outfp = fopen(outputFile, "w+"); 
     if(!outfp) { 
          fprintf(stderr, "Error opening \"%s\" for writing : %s\n", outputFile, strerror(errno));
          return false;
     }
     writePSPlotHeaderInfo(outfp);
     char *macroNames[4] = {"printdomainone", "printdomaintwo", "printdomainthree", "printdomainfour"}; 
     for(int i = 1; i <= 4; i++) { 
          char *domainDotBracketFile = runtimeConfig.getBranchTypeOutputFile((BranchID_t) i, "dot"); 
          fprintf(outfp, "%% Begin domain #%d structure printing routine: \n", i);
          fprintf(outfp, "/%s {\n\n", macroNames[i - 1]);
          fprintf(outfp, 
                  "   newpath\n"
                  "   /Helvetica findfont 32 scalefont setfont\n"
                  "   2 2 moveto\n"
                  "   setrgbcolor\n"
                  "   (Domain #%d: ) show\n\n", i);
          writePSPlotDomainData(outfp, domainDotBracketFile);
          fprintf(outfp, "} def\n%% End domain #%d structure printing routine --\n\n", i); 
          free(domainDotBracketFile);
     }
     writePSPlotActionData(outfp); 
     writePSPlotFooterInfo(outfp);
     fclose(outfp);  
     return true;
}

void Util::writePSPlotHeaderInfo(FILE *fp) { 
     fprintf(fp,
             "%%!PS-Adobe-3.0 EPSF-3.0\n"
             "%%%%Creator: RNABranchIDViz\n"
             "%%%%Title: RNA Secondary Structure Plot\n"
             "%%%%BoundingBox: 0 0 %d %d\n"
             "%%%%DocumentFonts: Helvetica\n"
             "%%%%Pages: 1\n"
             "%%%%EndComments\n\n", 2 * 700 + PSPLOT_DIVIDER_STROKE_SIZE, 2 * 700 + PSPLOT_DIVIDER_STROKE_SIZE);
     fprintf(fp, "%% to switch off outline pairs of sequence comment or\n"
             "%% delete the appropriate line near the end of the file\n\n");
     fprintf(fp, "%%%%BeginProlog\n");
     fprintf(fp, 
     "/RNAplot 100 dict def\n"
     "RNAplot begin\n"
     "/fsize  14 def\n"
     "%%/outlinecolor {0.2 setgray} bind def\n"
     "%%/paircolor    {0.2 setgray} bind def\n"
     "%%/seqcolor     {0   setgray} bind def'\n"
     "/outlinecolor {} bind def\n"
     "/paircolor    {} bind def\n"
     "/seqcolor     {} bind def\n"
     "/cshow  { dup stringwidth pop -2 div fsize -3 div rmoveto show} bind def\n"
     "/min { 2 copy gt { exch } if pop } bind def\n"
     "/max { 2 copy lt { exch } if pop } bind def\n"
     "/drawoutline {\n"
     "   gsave outlinecolor newpath\n"
     "   coor 0 get aload pop 0.8 0 360 arc %% draw 5' circle of 1st sequence\n"
     "   currentdict /cutpoint known        %% check if cutpoint is defined\n"
     "   {coor 0 cutpoint getinterval\n"
     "   {aload pop lineto} forall         %% draw outline of 1st sequence\n"
     "   coor cutpoint 1 add get aload pop\n"
     "   2 copy moveto 0.8 0 360 arc       %% draw 5' circle of 2nd sequence\n"
     "   coor cutpoint 1 add coor length cutpoint 1 add sub getinterval\n"
     "   {aload pop lineto} forall}        %% draw outline of 2nd sequence\n"
     "   {coor {aload pop lineto} forall}   %% draw outline as a whole\n"
     "   ifelse\n"
     "   stroke grestor\n"
     "} bind def\n"
     "/drawpairs {\n"
     "   paircolor\n"
     "   0.7 setlinewidth\n"
     "   [9 3.01] 9 setdash\n"
     "   newpath\n"
     "   pairs {aload pop\n"
     "       currentdict (cpr) known\n"
     "       { exch dup\n"
     "       coor  exch 1 sub get aload pop moveto\n"
     "       exch arccoords curveto\n"
     "       }\n"
     "       { coor exch 1 sub get aload pop moveto\n"
     "       coor exch 1 sub get aload pop lineto\n"
     "       } ifelse\n"
     "   } forall\n"
     "   stroke\n"
     "} bind def\n"
     "%% draw bases\n"
     "/drawbases {\n"
     "   [] 0 setdash\n"
     "   seqcolor\n"
     "   0\n"
     "   coor {\n"
     "       aload pop moveto\n"
     "       dup sequence exch 1 getinterval cshow\n"
     "       1 add\n"
     "   } forall\n"
     "   pop\n"
     "} bind def\n"
     "/init {\n"
     "   /Helvetica findfont fsize scalefont setfont\n"
     "   1 setlinejoin\n"
     "   1 setlinecap\n"
     "   0.8 setlinewidth\n"
     "   %% find the coordinate range\n"
     "   /xmax -1000 def /xmin 10000 def\n"
     "   /ymax -1000 def /ymin 10000 def\n"
     "   coor {\n"
     "      aload pop\n"
     "      dup ymin lt {dup /ymin exch def} if\n"
     "      dup ymax gt {/ymax exch def} {pop} ifelse\n"
     "      dup xmin lt {dup /xmin exch def} if\n"
     "      dup xmax gt {/xmax exch def} {pop} ifelse\n"
     "   } forall\n"
     "   /size {xmax xmin sub ymax ymin sub max} bind def\n"
     "   /width {xmax xmin sub} bind def\n"
     "   /height {ymax ymin sub} bind def\n"
     "   10 10 translate\n"
     "   680 size 10 add div dup scale\n"
     "   size width sub width xmin sub xmax sub add 2 div 5 add\n"
     "   size height sub height ymin sub ymax sub add 2 div 5 add\n"
     "   translate\n"
     "} bind def\n"
     "end\n");
     fprintf(fp, "%%%%EndProlog\n");
     fprintf(fp, "RNAplot begin\n"
             "%% data start here\n");
} 

// mostly modified from the ViennaRNA package's vrna_file_PS_rnaplot_a function:
void Util::writePSPlotDomainData(FILE *xyplot, const char *dotBracketSourceFile) { 

  FILE *fpDotFile = fopen(dotBracketSourceFile, "r"); 
  if(!fpDotFile) {
       fprintf(stderr, "Unable to open \"%s\" for reading : %s\n", dotBracketSourceFile, strerror(errno)); 
       return;
  }
  char *fileLine = NULL;
  size_t lineLength = 0;
  getline(&fileLine, &lineLength, fpDotFile); // gets the comment line 
  free(fileLine); 
  char *seq = NULL, *rec_rest = NULL;
  getline(&seq, &lineLength, fpDotFile); // gets the sequence line 
  if(seq[lineLength - 1] == '\n')
       seq[lineLength - 1] = '\0'; // kill the trailing newlinw
  getline(&rec_rest, &lineLength, fpDotFile); // gets the structure line 
  if(rec_rest[lineLength - 1] == '\n')
       rec_rest[lineLength - 1] = '\0';
  fclose(fpDotFile);

  char *structure = vrna_extract_record_rest_structure((const char **) &rec_rest, 0, 0); 
  vrna_md_t *md_p;
  vrna_md_set_default(md_p);

  float  xmin, xmax, ymin, ymax;
  int    i, length;
  int    ee, gb, ge, Lg, l[3];
  float *X, *Y;
  short *pair_table, *pair_table_g;
  char  *c, *string;
  vrna_md_t   md;

  if(!md_p){
    set_model_details(&md);
    md_p  = &md;
  }
    
  string = strdup(seq);
  length = strlen(string);

  pair_table = vrna_ptable(structure);
  pair_table_g = vrna_ptable(structure);

  ge=0;
  while ( (ee=parse_gquad(structure+ge, &Lg, l)) >0 ) {
    ge += ee;
    gb=ge-Lg*4-l[0]-l[1]-l[2]+1;
    /* add pseudo-base pair encloding gquad */
    for (i=0; i<Lg; i++) {
      pair_table_g[ge-i]=gb+i;
      pair_table_g[gb+i]=ge-i;
    }
  } 
      
  X = (float *) vrna_alloc((length+1)*sizeof(float));
  Y = (float *) vrna_alloc((length+1)*sizeof(float));
  i = naview_xy_coordinates(pair_table_g, X, Y);
  //if(i!=length)
  //  fprintf(stderr, "Warning: strange things happening in PS_rna_plot...");

  xmin = xmax = X[0];
  ymin = ymax = Y[0];
  for (i = 1; i < length; i++) {
     xmin = X[i] < xmin ? X[i] : xmin;
     xmax = X[i] > xmax ? X[i] : xmax;
     ymin = Y[i] < ymin ? Y[i] : ymin;
     ymax = Y[i] > ymax ? Y[i] : ymax;
  }

  /* cut_point */
  if ((c = strchr(structure, '&'))) {
    int cutpoint;
    cutpoint = c - structure;
    string[cutpoint] = ' '; /* replace & with space */
    fprintf(xyplot, "   /cutpoint %d def\n", cutpoint);
  }

  /* sequence */
  fprintf(xyplot,"   /sequence (\\\n");
  i=0;
  while (i<length) {
    fprintf(xyplot, "   %.255s\\\n", string+i);  /* no lines longer than 255 */
    i+=255;
  }
  fprintf(xyplot,"   ) def\n");
  /* coordinates */
  fprintf(xyplot, "   /coor [\n");
  for (i = 0; i < length; i++)
    fprintf(xyplot, "   [%3.8f %3.8f]\n", X[i], Y[i]);
  fprintf(xyplot, "   ] def\n");
  /* base pairs */
  fprintf(xyplot, "   /pairs [\n");
  for (i = 1; i <= length; i++)
    if (pair_table[i]>i)
      fprintf(xyplot, "   [%d %d]\n", i, pair_table[i]);
  /* add gquad pairs */
  ge=0;
  while ( (ee=parse_gquad(structure+ge, &Lg, l)) >0 ) {
    int k;
    fprintf(xyplot, "   %% gquad\n");
    ge += ee;
    gb=ge-Lg*4-l[0]-l[1]-l[2]+1; /* add pseudo-base pair encloding gquad */
    for (k=0; k<Lg; k++) {
      int ii, jj, il;
      for (il=0, ii=gb+k; il<3; il++) {
        jj = ii+l[il]+Lg;
        fprintf(xyplot, "   [%d %d]\n", ii, jj);
        ii = jj;
      }
      jj = gb+k;
      fprintf(xyplot, "   [%d %d]\n", jj, ii);
    }
  }

  fprintf(xyplot, "   ] def\n\n");
  fprintf(xyplot, "   init\n\n");
  /* draw the data */
  fprintf(xyplot,
          "   %% switch off outline pairs or bases by removing these lines\n"
          "   drawoutline\n"
          "   drawpairs\n"
          "   drawbases\n");

  free(seq);
  free(rec_rest); 
  free(structure);
  free(string);
  free(pair_table);
  free(pair_table_g);
  free(X); free(Y);
}

void Util::writePSPlotActionData(FILE *fp) { 
     fprintf(fp, "%% Now generate the output diagram:\n"); 
     fprintf(fp, 
             "%d setlinewidth\n"
             "0 700 newpath moveto\n"
             "%d 0 rlineto\n"
             "0 0 0 setrgbcolor\n"
             "stroke\n"
             "700 0 newpath moveto\n"
             "0 %d rlineto\n"
             "0 0 0 setrgbcolor\n"
             "stroke\n", PSPLOT_DIVIDER_STROKE_SIZE, 2 * 700 + PSPLOT_DIVIDER_STROKE_SIZE, 2 * 700 + PSPLOT_DIVIDER_STROKE_SIZE);
     fprintf(fp, 
             "initgraphics\n"
             "0 0 translate\n"
             "%0.3g %0.3g %0.3g printdomainone\n", 92.0 / 255, 160.0 / 255, 215.0 / 255);
     fprintf(fp, 
             "initgraphics\n"
             "0 %d translate\n"
             "%0.3g %0.3g %0.3g printdomaintwo\n", 700 + PSPLOT_DIVIDER_STROKE_SIZE, 183.0 / 255, 127.0 / 255, 77.0 / 255);
     fprintf(fp, 
             "initgraphics\n"
             "%d %d translate\n"
             "%0.3g %0.3g %0.3g printdomainthree\n", 700 + PSPLOT_DIVIDER_STROKE_SIZE, 700 + PSPLOT_DIVIDER_STROKE_SIZE, 243.0 / 255, 153.0 / 255, 193.0 / 255);
     fprintf(fp, 
             "initgraphics\n"
             "%d 0 translate\n"
             "%0.3g %0.3g %0.3g printdomainfour\n", 700 + PSPLOT_DIVIDER_STROKE_SIZE, 700 + PSPLOT_DIVIDER_STROKE_SIZE, 123.0 / 255, 204.0 / 255, 153.0 / 255);
     fprintf(fp, "\n"); 
} 

void Util::writePSPlotFooterInfo(FILE *fp) { 
     fprintf(fp, "%% Start Annotations\n");
     fprintf(fp, "%% End Annotations\n");
     fprintf(fp, "%% show it\nshowpage\n");
     fprintf(fp, "end\n");
     fprintf(fp, "%%%%EOF\n");
}

