/* DiagramImage.cpp 
   Author:  Maxie D. Schmidt (maxieds@gmail.com)
   Created: 2018.06.17
*/ 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "DiagramImage.h"
#include "BranchTypeIdentification.h"
#include "Utils.h"

#ifdef WITH_CAIRO_SUPPORT

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
     crSurface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, IMAGE_WIDTH, IMAGE_HEIGHT); 
     crDraw = cairo_create(crSurface);
}

DiagramImage_t::~DiagramImage_t() { 
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

    cairo_pattern_t *circleMask;
    cairo_push_group(crDraw);
    cairo_set_source_rgba (crDraw, 1.0, 1.0, 1.0, 0.0);
    cairo_paint(crDraw);
    cairo_arc(crDraw, IMAGE_WIDTH / 2, IMAGE_WIDTH / 2, radius, 0.0, 2.0 * M_PI);
    cairo_set_source_rgba (crDraw, 0.4, 0.4, 0.4, 1.0);
    cairo_fill (crDraw);
    circleMask = cairo_pop_group (crDraw);

    cairo_push_group(crDraw); 
    cairo_set_source_rgb(crDraw, 1.0, 1.0, 1.0); 
    cairo_paint(crDraw);
    cairo_set_line_width(crDraw, pixelWidth);

    for(unsigned int ui = 0; ui < numBases; ++ui) {
        const RNAStructure::BaseData* baseData1 = rnaStruct->GetBaseAt(ui);
	if (baseData1->m_pair != RNAStructure::UNPAIRED && baseData1->m_pair > ui) {
             RNABranchType_t::SetBranchColor(crDraw, rnaStruct->GetBranchTypeAt(ui)->getBranchID());
             DrawArc2(ui, baseData1->m_pair, centerX, centerY, angleBase, angleDelta, radius);
	     counter++;
        }
    }
    
    cairo_pop_group_to_source(crDraw);
    cairo_mask(crDraw, circleMask);
    cairo_pattern_destroy(circleMask);
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

void DiagramImage_t::DrawArc2(
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

    if (arc2 - arc1 > 180.0)
	arc1 += 360.0;
    if (arc1 - arc2 > 180.0)
	arc2 += 360.0;
    
    int boundingBoxCenterX = boundX + boundSize / 2;
    int boundingBoxCenterY = boundY + boundSize / 2;
    float boundingBoxRadius = boundSize / 2.0; //sqrt(2);

    if (arc2 > arc1) {
	cairo_arc(crDraw, boundingBoxCenterX, boundingBoxCenterY, boundingBoxRadius, arc1, arc2);
    }
    else {
	cairo_arc(crDraw, boundingBoxCenterX, boundingBoxCenterY, boundingBoxRadius, arc2, arc1);
    }
    cairo_stroke(crDraw); 

    cairo_move_to (crDraw, 0, 0);
    cairo_line_to (crDraw, 10, 10);
    cairo_stroke (crDraw);
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
    //angleDelta = (M_PI * 2.0f) / (float)numBases;
    //angleBase = 1.5f * M_PI;
    centerX = (float)resolution / 2.0f;
    centerY = (float)resolution / 2.0f;
    //centerX = (float)resolution / 4.0f;
    //centerY = (float)resolution / 4.0f;
    radius = centerX < centerY ? centerX - 15.f : centerY - 15.f;
}

#endif
