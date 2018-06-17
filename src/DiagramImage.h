/* DiagramImage.h : 
   Intermediate class used to create PNG images from the RNAStructure* data. 
   Based on RNAStructViz/DiagramWindow.*.
   Author:  Maxie D. Schmidt (maxieds@gmail.com)
   Created: 2018.06.17
*/ 

#ifndef __DIAGRAM_IMAGE_H__
#define __DIAGRAM_IMAGE_H__

#include <FL/Fl_Window.H>
#include <FL/x.H>

#include "RNAStructure.h"

#define IMAGE_WIDTH       (2048)
#define IMAGE_HEIGHT      (2048)
#define IMAGE_DEPTH       (3)

class DiagramImage_t : public Fl_Window { 

     private:
          RNAStructure *rnaStruct;
          int numPairs;
          unsigned char *pixelBuf;
          Fl_Offscreen offscreenImage;
          int pixelWidth;
          bool bufferDrawn;

     public:
          DiagramImage_t(RNAStructure *rnaStructBase); 
          ~DiagramImage_t();     

          int writePNGImage(const char *outputFilePath); 
     
     protected:
          void drawBuffers(); 

          void ComputeNumPairs(RNAStructure* structure);

          void ComputeCircle(
		const float& x1,
		const float& y1,
		const float& x2,
		const float& y2,
		const float& x3,
		const float& y3,
		double& cX,
		double& cY,
		double& r);

          void DrawArc(
		const unsigned int b1,
		const unsigned int b2,
		const float centerX,
		const float centerY,
		const float angleBase,
		const float angleDelta,
		const float radius);

         void DrawBase(
		const unsigned int index,
		const RNAStructure::Base base,
		const float centerX,
		const float centerY,
		const float angleBase,
		const float angleDelta,
		const float radius);

         void ComputeDiagramParams(
		const int numBases,
		const int resolution,
		float& centerX,
		float& centerY,
		float& angleBase,
		float& angleDelta,
		float& radius);

}; 

#endif
