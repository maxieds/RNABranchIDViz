/* DiagramImage.cpp 
   Author:  Maxie D. Schmidt (maxieds@gmail.com)
   Created: 2018.06.17
*/ 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <FL/fl_draw.H>
#include "DiagramImage.h"
#include "BranchTypeIdentification.h"
#include "Utils.h"

DiagramImage_t::DiagramImage_t(RNAStructure *rnaStructBase) : 
     rnaStruct(rnaStructBase), bufferDrawn(false) , numPairs(0) {
     if (rnaStruct->GetLength() > 1000) {
          pixelWidth = 1;
     } 
     else if(rnaStruct->GetLength() <= 1000 && rnaStruct->GetLength() >= 500) {
          pixelWidth = 2;
     } 
     else {
          pixelWidth = 3;
     }
     //pixelBuf = (unsigned char *) malloc(sizeof(unsigned char) * IMAGE_HEIGHT * IMAGE_WIDTH * IMAGE_DEPTH); 
     //memset(pixelBuf, 0xff, IMAGE_HEIGHT * IMAGE_WIDTH * IMAGE_DEPTH);
     crSurface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, IMAGE_WIDTH, IMAGE_HEIGHT); 
     crDraw = cairo_create(crSurface);
}

DiagramImage_t::~DiagramImage_t() { 
     //free(pixelBuf);
     //pixelBuf = NULL;
     cairo_destroy(crDraw);
     cairo_surface_destroy(crSurface);
} 

int DiagramImage_t::writePNGImage(const char *outputFilePath) { 
     if(!bufferDrawn)
          drawBuffers(); 
     return cairo_surface_write_to_png(crSurface, outputFilePath); 
}

void DiagramImage_t::drawBuffers() { 
    
    float centerX = 0.0f;
    float centerY = 0.0f;
    float angleBase = 0.0f;
    float angleDelta = 0.0f;
    float radius = 0.0f;
    int counter = 0;
    unsigned int numBases = rnaStruct->GetLength();
    
    ComputeNumPairs(rnaStruct); 
    ComputeDiagramParams(numBases, IMAGE_WIDTH, centerX, centerY, angleBase, angleDelta, radius);

    cairo_set_source_rgb(crDraw, 1.0, 1.0, 1.0); 
    cairo_paint(crDraw);
    cairo_set_line_width(crDraw, pixelWidth);

    for(unsigned int ui = 0; ui < numBases; ++ui) {
        const RNAStructure::BaseData* baseData1 = rnaStruct->GetBaseAt(ui);
	if (baseData1->m_pair != RNAStructure::UNPAIRED && baseData1->m_pair > ui) {
             RNABranchType_t::SetBranchColor(crDraw, rnaStruct->GetBranchTypeAt(ui).getBranchID());
             DrawArc(ui, baseData1->m_pair, centerX, centerY, angleBase, angleDelta, radius);
	     counter++;
        }
    }
    bufferDrawn = true;

}

void DiagramImage_t::ComputeNumPairs(RNAStructure *structure) {
    unsigned int numBases = structure->GetLength();
    unsigned int counter1 = 0;
    for(unsigned int ui = 0; ui < numBases; ++ui) {
         const RNAStructure::BaseData* baseData1 = rnaStruct->GetBaseAt(ui);
         if(baseData1->m_pair != RNAStructure::UNPAIRED && baseData1->m_pair > ui)
              counter1++;
    }
    numPairs = counter1;
}

void DiagramImage_t::ComputeCircle(
    const float& x1,
    const float& y1,
    const float& x2,
    const float& y2,
    const float& x3,
    const float& y3,
    double& cX,
    double& cY,
    double& r)
{
    double denom = x1 * (y2 - y3) - y1 * (x2 - x3) + x2 * y3 - y2 * x3;

    if (denom < 0.001)
    {
		cX = cY = 0.0f;
		r = 0.0f;
    }

    double sq1 = x1 * x1 + y1 * y1;
    double sq2 = x2 * x2 + y2 * y2;
    double sq3 = x3 * x3 + y3 * y3;

    cX = (sq1 * (y2 - y3) - y1 * (sq2 - sq3) + sq2 * y3 - y2 * sq3) / (2.0 * 
    	denom);
    cY = (x1 * (sq2 - sq3) - sq1 * (x2 - x3) + sq3 * x2 - x3 * sq2) / (2.0 * 
    	denom);

    r = sqrt((x1 - cX) * (x1 - cX) + (y1 - cY) * (y1 - cY));
}

void DiagramImage_t::DrawArc(
    const unsigned int b1,
    const unsigned int b2,
    const float centerX,
    const float centerY,
    const float angleBase,
    const float angleDelta,
    const float radius)
{
    
    float angle1 = angleBase - (float)b1 * angleDelta;
    float xPosn1 = centerX + cos(angle1) * radius;
    float yPosn1 = centerY - sin(angle1) * radius;

    float angle2 = angleBase - (float)b2 * angleDelta;
    float xPosn2 = centerX + cos(angle2) * radius;
    float yPosn2 = centerY - sin(angle2) * radius;

    // Calculate a third point on the arc, midway between the endpoints.
    float midAngle = (angle1 + angle2) / 2.0f;
    float diffAngleRatio = (angle1 - angle2) / M_PI;
    float xPosn3 = centerX + cos(midAngle) * radius * (1.0f - diffAngleRatio);
    float yPosn3 = centerY - sin(midAngle) * radius * (1.0f - diffAngleRatio);

    double arcX = 0.0f;
    double arcY = 0.0f;
    double arcR = 0.0f;
    ComputeCircle(xPosn1, yPosn1, xPosn2, yPosn2, xPosn3, yPosn3, arcX, arcY, arcR);

    int boundX = (int)(arcX - arcR);
    int boundY = (int)(arcY - arcR);
    int boundSize = (int)(2.0f * arcR);
    double arc1 = 180.0 / M_PI * atan2(arcY - yPosn1, xPosn1 - arcX);
    double arc2 = 180.0 / M_PI * atan2(arcY - yPosn2, xPosn2 - arcX);

    //boundSize = boundSize * 0.7071;
    //cairo_translate(crDraw, boundX, boundY); 
    //cairo_scale(crDraw, boundSize, boundSize);
    if (arc2 - arc1 > 180.0)
	arc1 += 360.0;
    if (arc1 - arc2 > 180.0)
	arc2 += 360.0;
    if (arc2 > arc1) {
	cairo_arc(crDraw, boundX, boundY, boundSize, arc1, arc2);
        //cairo_arc(crDraw, boundX, boundY, boundSize, arc1, arc2);
        //if (pixelWidth > 1) 
        //    cairo_arc(crDraw, boundX + 1, boundY + 1, boundSize, arc1, arc2);
        //if (pixelWidth > 2)
        //    cairo_arc(crDraw, boundX + 1, boundY, boundSize, arc1, arc2);
    }
    else {
	cairo_arc(crDraw, boundX, boundY, boundSize, arc2, arc1);
        //cairo_arc(crDraw, boundX, boundY, boundSize, arc2, arc1);
        //if (pixelWidth > 1)
        //    cairo_arc(crDraw, boundX + 1, boundY + 1, boundSize, arc2, arc1);
        //if (pixelWidth > 2)
        //    cairo_arc(crDraw, boundX + 1, boundY, boundSize, arc2, arc1);
    }
    cairo_stroke(crDraw); 

}

void DiagramImage_t::DrawArc2(
    const unsigned int b1,
    const unsigned int b2,
    const float centerX,
    const float centerY,
    const float angleBase,
    const float angleDelta,
    const float radius)
{
    
    int pairDistance = MAX(b1, b2) - MIN(b1, b2); 
    float h = radius * pairDistance / numPairs;
    float theta1 = angleBase - angleDelta * MIN(b1, b2);
    float theta2 = angleBase - angleDelta * MAX(b1, b2);
    float P1x = radius * cos(theta1), P1y = radius * sin(theta1); 
    float P2x = radius * cos(theta2), P2y = radius * sin(theta2); 
    float arc1 = -1 * acos(P1x / h), arc2 = asin(P2y / h); 
    float x0, y0;
    x0 = ((P1x + P2x)*(pow(P1x - P2x,2) + pow(P1y - P2y,2)) - sqrt((-pow(P1x - P2x,2) - pow(P1y - P2y,2)) *
          (-4*pow(h,2) + pow(P1x - P2x,2) + pow(P1y - P2y,2))*pow(P1y - P2y,2)))/(2.0*(pow(P1x - P2x,2) + pow(P1y - P2y,2)));
    y0 = (pow(P1y,4) + pow(P1y,2)*pow(P1x - P2x,2) + P1x*sqrt((-pow(P1x - P2x,2) - pow(P1y - P2y,2))*
         (-4*pow(h,2) + pow(P1x - P2x,2) + pow(P1y - P2y,2))*pow(P1y - P2y,2))- P2x*sqrt((-pow(P1x - P2x,2) - pow(P1y - P2y,2))*
         (-4*pow(h,2) + pow(P1x - P2x,2) + pow(P1y - P2y,2))*pow(P1y - P2y,2))- 2*pow(P1y,3)*P2y + 2*P1y*pow(P2y,3) - 
         pow(P2y,2)*(pow(P1x - P2x,2) + pow(P2y,2)))/ (2.0*(pow(P1x - P2x,2) + pow(P1y - P2y,2))*(P1y - P2y));
    if (arc2 - arc1 > 180.0)
	arc1 += 360.0;
    if (arc1 - arc2 > 180.0)
	arc2 += 360.0;

    cairo_arc(crDraw, x0, y0, h, MIN(arc1, arc2), MAX(arc1, arc2));
    cairo_stroke(crDraw); 
}

void DiagramImage_t::ComputeDiagramParams(
    const int numBases,
    const int resolution,
    float& centerX,
    float& centerY,
    float& angleBase,
    float& angleDelta,
    float& radius)
{
    angleDelta = (M_PI * 2.0f - 0.05f) / (float)numBases;
    angleBase = 1.5f * M_PI - 0.025f;
    centerX = (float)resolution / 2.0f;
    centerY = (float)resolution / 2.0f;
    radius = centerX < centerY ? centerX - 15.f : centerY - 15.f;
}

