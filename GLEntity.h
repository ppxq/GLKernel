/////////////////////////////////////////////////////////////////////////////
//
//	GLEntity - the basic class of entity can be displayed by GLKernel
//
//
/////////////////////////////////////////////////////////////////////////////


#ifndef CW_GLENTITY
#define CW_GLENTITY

#include <afx.h>
#include <gl\GL.h>
#include "GL/glu.h"
#include "GL/glaux.h"

class GLEntity : public CObject  
{
public:
	GLEntity() {bShow=TRUE; entityType=0;};
	virtual ~GLEntity() {};

	// Implement the virtual method which draw this entity
	//		TRUE - draw the shading mode
	//		FALSE - draw the mesh mode
	virtual void drawShade() {};
	virtual void drawProfile() {};
	virtual void drawPreMesh() {};
	virtual void drawMesh() {};
	virtual void drawHighLight() {};
	virtual void drawPick() {};
	virtual void drawNode() {};
	virtual void drawNodeNormal(const double &length) {};
	virtual void drawLineStrip() {};
	virtual void drawEdgeNode() {};

	// Implement the maximum distance to the original point of this entity 
	virtual float getRange() {return 0.0;};

	BOOL bShow;

	UINT entityType;

protected:
	float red, green, blue;
};

#endif
