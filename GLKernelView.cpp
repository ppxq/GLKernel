/////////////////////////////////////////////////////////////////////////////
//	CGLKernelView - the kernel class of OpenGL display & manipulation
//
/////////////////////////////////////////////////////////////////////////////

// GLKernelView.cpp : implementation file
//

#include "stdafx.h"
#include "GLKernelView.h"
#include "GL/glut.h"
#include <GL/gl.h>
#include <GL/glaux.h>
#include <GL/glu.h>
#include <math.h>

#include "MouseTool.h"
#include "GLList.h"
#include "CW_Geometry.h"
#include "Position_array.h"



#pragma warning(disable : 4244)     // MIPS
#pragma warning(disable : 4136)     // X86
#pragma warning(disable : 4051)     // ALPHA

#define DEGREE_TO_RATATE(x)		0.0174532922222*x
#define BUFSIZE					2048

/////////////////////////////////////////////////////////////////////////////
// CGLKernelView
CGLKernelView::CGLKernelView(CWnd *pWnd)
{
	m_pWnd = pWnd;
	m_currentTool=NULL;
	m_displayObjList.RemoveAll();
	m_HighLightObj=NULL;
	m_Shading=TRUE;
	m_Mesh=FALSE;
	m_Profile=FALSE;
	m_Node=FALSE;
	m_NodeNormal=FALSE;

	ShadingType=1;

}

CGLKernelView::~CGLKernelView()
{
	clear_tools();
	clear_gllist();

	POSITION Pos;
	for(Pos=m_displayObjList.GetHeadPosition();Pos!=NULL;)
	{
		GLEntity *entity=(GLEntity *)(m_displayObjList.GetNext(Pos));
		if (entity) {delete entity;	entity=NULL;}
	}
	ClearDisplayObjList();


}

/* -------------------------------- Public Functions -------------------------------*/
void CGLKernelView::ClearAll()
{
	clear_tools();
	clear_gllist();

	POSITION Pos;
	for(Pos=m_displayObjList.GetHeadPosition();Pos!=NULL;)
	{
		GLEntity *entity=(GLEntity *)(m_displayObjList.GetNext(Pos));
		if (entity) delete entity;
	}
	ClearDisplayObjList();

	InitValue();

	refresh();
}

void CGLKernelView::OnDraw(CDC* pDC)
{
	// TODO: add draw code here
	refresh();
}

BOOL CGLKernelView::PreCreateWindow(CREATESTRUCT& cs)
{
	cs.style |= WS_CLIPCHILDREN | WS_CLIPSIBLINGS;

	return CView::PreCreateWindow(cs);
}

int CGLKernelView::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{
	// TODO: Add your specialized creation code here
	PIXELFORMATDESCRIPTOR pfd =
    {
        sizeof(PIXELFORMATDESCRIPTOR), // Structure size.
        1,                             // Structure version number.
        PFD_DRAW_TO_WINDOW |           // Property flags.
		PFD_DOUBLEBUFFER   |
        PFD_SUPPORT_OPENGL,
        PFD_TYPE_RGBA,
        24,                            // 24-bit color.
        0, 0, 0, 0, 0, 0,              // Not concerned with these.
        0, 0, 0, 0, 0, 0, 0,           // No alpha or accum buffer.
        32,                            // 32-bit depth buffer.
        0, 0,                          // No stencil or aux buffer.
        PFD_MAIN_PLANE,                // Main layer type.
        0,                             // Reserved.
        0, 0, 0                        // Unsupported.
    };


	CDC *pDC = m_pWnd->GetDC();

    int pixelFormat = ChoosePixelFormat(pDC->m_hDC, &pfd);
    SetPixelFormat(pDC->m_hDC, pixelFormat, &pfd);
	DescribePixelFormat(pDC->m_hDC, pixelFormat, sizeof(pfd), &pfd);
	m_hRC = wglCreateContext(pDC->m_hDC);

	InitValue();
	return 0;
}

void CGLKernelView::OnDestroy() 
{
	CView::OnDestroy();

	wglMakeCurrent(NULL,NULL);
	wglDeleteContext(m_hRC);
}

void CGLKernelView::OnSize(UINT nType, int cx, int cy) 
{
	CView::OnSize(nType, cx, cy);
	
	// TODO: Add your message handler code here
	CDC *pDC = m_pWnd->GetWindowDC();
	if (!pDC) return;

	m_SizeX=cx;		m_SizeY=cy;

	wglMakeCurrent(pDC->m_hDC,m_hRC);

	GLInit();

    wglMakeCurrent(NULL,NULL);

	m_pWnd->ReleaseDC(pDC);

	refresh();
}

void CGLKernelView::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags) 
{
	// TODO: Add your message handler code here and/or call default
	pick_event pe;
	pe.nChar=nChar;
	pe.nFlags=nFlags;

	if (nChar==VK_CANCEL) 
		MessageBox("This GLKernal is developed by Charlie C. L. WANG!", "Copyright");

	if (m_currentTool) m_currentTool->process_event(KEY_PRESS,pe);
	
	CView::OnChar(nChar, nRepCnt, nFlags);
}

void CGLKernelView::OnMouseMove(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	pick_event pe;
	pe.x=(double)point.x;
	pe.y=(double)point.y;
	pe.nFlags=nFlags;

	if (m_currentTool) m_currentTool->process_event(MOUSE_MOVE,pe);
	
	CView::OnMouseMove(nFlags, point);
}

void CGLKernelView::OnLButtonUp(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	pick_event pe;
	pe.x=(double)point.x;
	pe.y=(double)point.y;
	pe.nFlags=nFlags;

	if (m_currentTool) m_currentTool->process_event(MOUSE_BUTTON_UP,pe);
	
	CView::OnLButtonUp(nFlags, point);
}

void CGLKernelView::OnLButtonDown(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	pick_event pe;
	pe.x=(double)point.x;
	pe.y=(double)point.y;
	pe.nFlags=nFlags;

	if (m_currentTool) m_currentTool->process_event(MOUSE_BUTTON_DOWN,pe);
	//clear_tools();
	
	
	CView::OnLButtonDown(nFlags, point);
}

void CGLKernelView::OnRButtonUp(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	pick_event pe;
	pe.x=(double)point.x;
	pe.y=(double)point.y;
	pe.nFlags=nFlags;

	if (m_currentTool) m_currentTool->process_event(MOUSE_BUTTON_UP,pe);
	
	CView::OnRButtonUp(nFlags, point);
}

void CGLKernelView::OnRButtonDown(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	pick_event pe;
	pe.x=(double)point.x;
	pe.y=(double)point.y;
	pe.nFlags=nFlags;

	if (m_currentTool) m_currentTool->process_event(MOUSE_BUTTON_DOWN,pe);
	
	CView::OnRButtonDown(nFlags, point);
}

/* -------------------------------- Private Functions ------------------------------*/
void CGLKernelView::InitValue()
{
	m_xRotation = 0.0f;
	m_yRotation = 0.0f;

	m_xTranslation = 0.0f;
	m_yTranslation = 0.0f;
	m_zTranslation = 0.0f;

	m_Scaling = 0.5f;
	m_Range = 1.0f;

	//background color
	m_ClearColorRed   = 1.0;
	m_ClearColorGreen = 1.0;
	m_ClearColorBlue  = 1.0;
	/*m_ClearColorRed   = 0.35f;
	m_ClearColorGreen = 0.35f;
	m_ClearColorBlue  = 0.35f;*/

	m_red = 0.0;
	m_green = 1.0;
	m_blue = 0.0;
	
	m_lineWidth = 1;

	m_axisDisplay = true;
	m_lightON = true;
	RedrawAll();
}

void CGLKernelView::RedrawAll()
{
	m_shadeDrawn = false;
	m_nodeDrawn = false;
	m_nodeNormalDrawn = false;
}

void CGLKernelView::GLInit()
{
	int cx=m_SizeX;
	int cy=m_SizeY;
	float scale=m_Scaling*2.0f;

	if ((m_Range*scale)<0.5) scale=0.5/m_Range;

	glViewport(0,0,cx,cy);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	
    if (cx <= cy)
	{
	    glOrtho (-m_Range, m_Range, -m_Range*(GLfloat)cy/(GLfloat)cx, 
			m_Range*(GLfloat)cy/(GLfloat)cx, 
			-m_Range*scale, m_Range*scale);
		m_MappingScale=cx/(m_Range*2.0);
	}
    else 
	{
		glOrtho (-m_Range*(GLfloat)cx/(GLfloat)cy, 
			m_Range*(GLfloat)cx/(GLfloat)cy, -m_Range, m_Range, 
			-m_Range*scale, m_Range*scale);
		m_MappingScale=cy/(m_Range*2.0);
	}

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void CGLKernelView::GLEnableSilverShading(){
	glShadeModel(GL_SMOOTH);
	glEnable(GL_NORMALIZE);
	
	glEnable(GL_DEPTH_TEST);
    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL); 


	GLfloat	ambientProperties[] = {0.7f, 0.7f, 0.7f, 1.0f};
	GLfloat	diffuseProperties[]  = {0.8f, 0.8f, 0.8f, 1.0f};
	GLfloat	specularProperties[] = {1.0f, 1.0f, 1.0f, 1.0f};

	glLightfv( GL_LIGHT0, GL_AMBIENT, ambientProperties);
	glLightfv( GL_LIGHT0, GL_DIFFUSE, diffuseProperties);
	glLightfv( GL_LIGHT0, GL_SPECULAR, specularProperties);


	//silver	
	GLfloat ambientMaterial[]={0.19225, 0.19225, 0.19225};
	GLfloat diffuseMaterial[]={0.50754, 0.50754, 0.50754};
	GLfloat specularMaterial[]={0.508273, 0.508273, 0.508273};
	GLfloat shininessMaterial = 0.4;

	glMaterialfv(GL_FRONT, GL_AMBIENT, ambientMaterial);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuseMaterial);
	glMaterialfv(GL_FRONT, GL_SPECULAR, specularMaterial);
	glMaterialf(GL_FRONT, GL_SHININESS, shininessMaterial * 128.0);

	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);

	glEnable(GL_LIGHT0);
	//glEnable(GL_LIGHT1);
	glEnable(GL_LIGHTING);
}

void CGLKernelView::GLEnableGrayShading(){ //default opengl shading
	// Lights, material properties
	glShadeModel(GL_SMOOTH);
	glEnable(GL_NORMALIZE);
	
	glEnable(GL_DEPTH_TEST);
    glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL); 

	GLfloat	ambientProperties[]  = {0.0f, 0.0f, 0.0f, 1.0f};
	GLfloat	diffuseProperties[]  = {1.0f, 1.0f, 1.0f, 1.0f};
	GLfloat	specularProperties[] = {1.0f, 1.0f, 1.0f, 1.0f};

	glLightfv( GL_LIGHT0, GL_AMBIENT, ambientProperties);
	glLightfv( GL_LIGHT0, GL_DIFFUSE, diffuseProperties);
	glLightfv( GL_LIGHT0, GL_SPECULAR, specularProperties);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);

	GLfloat ambientMaterial[]={0.2, 0.2, 0.2};
	GLfloat diffuseMaterial[]={0.8, 0.8, 0.8};
	GLfloat specularMaterial[]={0.0, 0.0, 0.0};
	GLfloat shininessMaterial = 0.0;

	glMaterialfv(GL_FRONT, GL_AMBIENT, ambientMaterial);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuseMaterial);
	glMaterialfv(GL_FRONT, GL_SPECULAR, specularMaterial);
	glMaterialf(GL_FRONT, GL_SHININESS, shininessMaterial * 128.0);

	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
}

void CGLKernelView::GLEnableLight()
{
	switch(ShadingType){
		case 0: default:
			GLEnableGrayShading();
			//printf("gray shading\r");
			break;
		case 1:
			GLEnableSilverShading(); 
			//printf("silver shading\r");
			break;
	}
}

void CGLKernelView::GLDisableLight()
{
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHTING);
//	glDisable(GL_DEPTH_TEST);
}

void CGLKernelView::GLDrawAxis()
{
	double scale = m_Range/m_Scaling;
	double axisLength=0.1*scale;
	GLUquadricObj *quadratic = gluNewQuadric();     // Create A Pointer The Quadric Object ( NEW )

	gluQuadricNormals(quadratic, GLU_SMOOTH);       // Create Smooth Normals ( NEW )

	GLEnableLight();
	glPushMatrix();
	glTranslated(-1.25*m_Range,-0.85*m_Range,m_Range*m_Scaling);
	glRotatef(m_xRotation,1.0f,0.0f,0.0f);
	glRotatef(m_yRotation,0.0f,1.0f,0.0f);
	glScalef(m_Scaling,m_Scaling,m_Scaling);
	//draw axis
	glColor3f(1.0,0.0,0.0);         //      x-axis
	glPushMatrix();
	glRotatef(90.0, 0.0, 1.0, 0.0);
	gluCylinder(quadratic,0.01f*scale,0.01f*scale,axisLength,8,1);
	glTranslated(0.0, 0.0, axisLength);
	gluCylinder(quadratic,0.018f*scale,0.0f,0.02f*scale,8,1);       // Draw Our Cylinder
	glPopMatrix();
	glColor3f(0.0,1.0,0.0);         //      y-axis
	glPushMatrix();
	glRotatef(-90.0, 1.0, 0.0, 0.0);
	gluCylinder(quadratic,0.01f*scale,0.01f*scale,axisLength,8,1);
	glTranslated(0.0, 0.0, axisLength);
	gluCylinder(quadratic,0.018f*scale,0.0f,0.02f*scale,8,1);       // Draw Our Cylinder
	glPopMatrix();
	glColor3f(0.0,0.0,1.0);         //      z-axis
	glPushMatrix();
	gluCylinder(quadratic,0.01f*scale,0.01f*scale,axisLength,8,1);
	glTranslated(0.0, 0.0, axisLength);
	gluCylinder(quadratic,0.018f*scale,0.0f,0.02f*scale,8,1);       // Draw Our Cylinder
	glPopMatrix();
	glPopMatrix();
	GLDisableLight();
	gluDeleteQuadric(quadratic);                            // Delete Quadratic - Free Resources
}

void CGLKernelView::GLDrawPolyCount()
{
}

void CGLKernelView::GLDraw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(m_ClearColorRed,m_ClearColorGreen,m_ClearColorBlue,1.0f);

	////////////////////////////////////////////////////////////////
	//	The following lines are drawing axis.
	if (m_axisDisplay) GLDrawAxis();

	// Position / translation / scale
	glTranslated(m_xTranslation,m_yTranslation,m_zTranslation);
	glRotatef(m_xRotation,1.0f,0.0f,0.0f);
	glRotatef(m_yRotation,0.0f,1.0f,0.0f);
	glScalef(m_Scaling,m_Scaling,m_Scaling);

	//GLDisableLight();
	

	////////////////////////////////////////////////////////////////
	//	The following lines are drawing poly count.
	if (m_polyCountDisplay) GLDrawPolyCount();

	////////////////////////////////////////////////////////////////
	//	The following lines are drawing Object.
	//		Default rendering 
	glColor3f(0.7f,0.7f,0.7f);
	GLDrawDisplayObjList();

	////////////////////////////////////////////////////////////////
	//	The following lines are drawing GLList and MouseTools.
	GLDrawGLList();
	if (m_currentTool) m_currentTool->draw();

	glGetDoublev(GL_MODELVIEW_MATRIX, modelMatrix);
	glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
	glGetIntegerv(GL_VIEWPORT, viewport);
}

void CGLKernelView::GLDrawGLList()
{
	POSITION Pos;
	for(Pos=m_glList.GetHeadPosition();Pos!=NULL;)
	{
		GLList *glList=(GLList *)(m_glList.GetNext(Pos));
		if (!glList) continue;
		glList->draw(this);
	}
}

void CGLKernelView::GLDrawDisplayObjList()
{
	glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(0.5,0.5);
	glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
	glEnable(GL_POLYGON_OFFSET_LINE);
	glPolygonOffset(1.0,1.0);
	

	//GLDisableLight();
	GLEnableLight();

	POSITION Pos;

	if (m_Node)
	{
		if (!m_nodeDrawn)
		{
			glDeleteLists(2, 1);
			glNewList(2,GL_COMPILE_AND_EXECUTE);
			for(Pos=m_displayObjList.GetHeadPosition();Pos!=NULL;)
			{
				GLEntity *entity=(GLEntity *)(m_displayObjList.GetNext(Pos));
				if (!entity) continue;
				if (!(entity->bShow)) continue;
				entity->drawNode();
				m_nodeDrawn = true;
			}
			glEndList();
		} else glCallList(2);
	}

	if (m_lightON)
		GLEnableLight();
	else 
		GLDisableLight();

	if (m_Shading)
	{
		if (!m_shadeDrawn)
		{
			glDeleteLists(1, 1);
			glNewList(1,GL_COMPILE_AND_EXECUTE);
			for(Pos=m_displayObjList.GetHeadPosition();Pos!=NULL;)
			{
				GLEntity *entity=(GLEntity *)(m_displayObjList.GetNext(Pos));
				if (!entity) continue;
				if (!(entity->bShow)) continue;
				//if (entity==m_HighLightObj) 
				//	entity->drawHighLight();
				//else
					entity->drawShade();
				m_shadeDrawn = true;
			}
			glEndList();
		} else glCallList(1);
	}
	else
	{
		if (m_Mesh)
		{
			GLDisableLight();
			glEnable(GL_DEPTH_TEST);
			for(Pos=m_displayObjList.GetHeadPosition();Pos!=NULL;)
			{
				GLEntity *entity=(GLEntity *)(m_displayObjList.GetNext(Pos));
				if (!entity) continue;
				if (!(entity->bShow)) continue;
				entity->drawPreMesh();
			}
		}
	}

	GLDisableLight();
	if (m_NodeNormal)
	{
		if (!m_nodeNormalDrawn)
		{
			glDeleteLists(3, 1);
			glNewList(3,GL_COMPILE_AND_EXECUTE);
			for(Pos=m_displayObjList.GetHeadPosition();Pos!=NULL;)
			{
				GLEntity *entity=(GLEntity *)(m_displayObjList.GetNext(Pos));
				if (!entity) continue;
				if (!(entity->bShow)) continue;
				entity->drawNodeNormal(0.05*m_Range);
				m_nodeNormalDrawn = true;
			}
			glEndList();
		} else glCallList(3);
	}

	if (m_Mesh)
	{
		for(Pos=m_displayObjList.GetHeadPosition();Pos!=NULL;)
		{
			GLEntity *entity=(GLEntity *)(m_displayObjList.GetNext(Pos));
			if (!entity) continue;
			if (!(entity->bShow)) continue;
			/*if (entity==m_HighLightObj) 
				entity->drawHighLight();
			else*/
				entity->drawMesh();
		}
	}

	if (m_Profile)
	{
		for(Pos=m_displayObjList.GetHeadPosition();Pos!=NULL;)
		{
			GLEntity *entity=(GLEntity *)(m_displayObjList.GetNext(Pos));
			if (!entity) continue;
			if (!(entity->bShow)) continue;
			entity->drawProfile();
			/*if (entity==m_HighLightObj) 
				entity->drawHighLight();*/
		}
	}

}

void CGLKernelView::GLPickInit()
{
	int cx=m_SizeX;
	int cy=m_SizeY;
	float scale=m_Scaling;

	if ((m_Range*scale)<0.5) scale=0.5/m_Range;
    if (cx <= cy)
	{
	    glOrtho (-m_Range, m_Range, -m_Range*(GLfloat)cy/(GLfloat)cx, 
			m_Range*(GLfloat)cy/(GLfloat)cx, -m_Range*scale, m_Range*scale);
		m_MappingScale=cx/(m_Range*2.0);
	}
    else 
	{
		glOrtho (-m_Range*(GLfloat)cx/(GLfloat)cy, 
			m_Range*(GLfloat)cx/(GLfloat)cy, -m_Range, m_Range, -m_Range*scale, m_Range*scale);
		m_MappingScale=cy/(m_Range*2.0);
	}
}

void CGLKernelView::GLPickDraw(UINT nType)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(m_ClearColorRed,m_ClearColorGreen,m_ClearColorBlue,1.0f);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	glInitNames();
	glPushName(0);

	glTranslated(m_xTranslation,m_yTranslation,m_zTranslation);
	double ca=cos(DEGREE_TO_RATATE(m_yRotation));
	double sa=sin(DEGREE_TO_RATATE(m_yRotation));
	double cb=cos(DEGREE_TO_RATATE(m_xRotation));
	double sb=sin(DEGREE_TO_RATATE(m_xRotation));
	GLdouble R[16]={
		ca,sa*sb,-sa*cb,0,
		0,cb,sb,0,
		sa,-ca*sb,ca*cb,0,
		0,0,0,1
	};
	glMultMatrixd(R);
	glScalef(m_Scaling,m_Scaling,m_Scaling);

	GLDisableLight();
	////////////////////////////////////////////////////////////////
	//	The following lines are drawing Object.
	//		Default rendering 
	glColor3f(0.7f,0.7f,0.7f);
	GLPickDrawDisplayObjList(nType);

	glPopMatrix();
	glFlush();
}

void CGLKernelView::GLPickDrawDisplayObjList(UINT nType)
{
	UINT nameIndex=1;
	POSITION Pos;
	for(Pos=m_displayObjList.GetHeadPosition();Pos!=NULL;)
	{
		GLEntity *entity=(GLEntity *)(m_displayObjList.GetNext(Pos));

		if (!entity) continue;
		if ((entity->bShow) && (entity->entityType==nType))
		{
			glPushMatrix();
			glLoadName(nameIndex++);
			if (m_Shading)
			{
				GLEnableLight();
				entity->drawShade();
				GLDisableLight();
			}
			if (m_Mesh)
				entity->drawMesh();
			if (m_Profile)
				entity->drawProfile();
			glPopMatrix();
		}
	}
}

GLEntity* CGLKernelView::GetPickObjByName(UINT nType,UINT name)
{
	UINT nameIndex=1;

	POSITION Pos;
	for(Pos=m_displayObjList.GetHeadPosition();Pos!=NULL;)
	{
		GLEntity *entity=(GLEntity *)(m_displayObjList.GetNext(Pos));

		if (!entity) continue;
		if ((entity->bShow) && (entity->entityType==nType))
		{
			if (name==nameIndex) return entity;
			nameIndex++;
		}
	}

	return NULL;
}


//////////////////////////////////////////////////////////////////////////////////////
//
//
/* -------------------------------- Command Functions ------------------------------*/
//
//
//////////////////////////////////////////////////////////////////////////////////////
void CGLKernelView::GetScale(float &scale) 
{
	scale=m_Scaling;
}

void CGLKernelView::SetScale(float scale) 
{
	m_Scaling=scale;
}

void CGLKernelView::SetViewDirection(UINT nDirID)
{
	switch(nDirID)
	{
	case 0:{	//VD_FRONTVIEW
				m_xRotation=0.0;	m_yRotation=0.0;
				refresh();
		   }break;
	case 1:{	//VD_LEFTVIEW
				m_xRotation=0.0;	m_yRotation=90.0;
				refresh();
		   }break;
	case 2:{	//VD_RIGHTVIEW
				m_xRotation=0.0;	m_yRotation=-90.0;
				refresh();
		   }break;
	case 3:{	//VD_BACKVIEW
				m_xRotation=0.0;	m_yRotation=180.0;
				refresh();
		   }break;
	case 4:{	//VD_TOPVIEW	
				m_xRotation=90.0;	m_yRotation=0.0;
				refresh();
		   }break;
	case 5:{	//VD_BOTTOMVIEW
				m_xRotation=-90.0;	m_yRotation=0.0;
				refresh();
		   }break;
	case 6:{	//VD_ISOMETRICVIEW
				m_xRotation=27.0;	m_yRotation=-45.0;
				refresh();
		   }break;
	case 7:{	//VD_BACKISOMETRICVIEW
				m_xRotation=27.0;	m_yRotation=135.0;
				refresh();
		   }break;
	}
}


void CGLKernelView::refresh()
{
	CDC *pDC1 = m_pWnd->GetWindowDC();
	if (!pDC1) return;
	wglMakeCurrent(pDC1->m_hDC,m_hRC);
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glPushMatrix();
	GLInit();
	GLDraw();
	glPopMatrix();
	glFlush();
	SwapBuffers(pDC1->m_hDC);
	wglMakeCurrent(NULL, NULL);
	
}

void CGLKernelView::zoom_all_in_view()
{
	float newRange,oldRange=m_Range;
	m_Range=1.0f;
//
	POSITION Pos;	BOOL flag=TRUE;
	for(Pos=m_displayObjList.GetHeadPosition();Pos!=NULL;)
	{
		GLEntity *tempEntity=(GLEntity *)(m_displayObjList.GetNext(Pos));
		newRange=tempEntity->getRange();
		if ((newRange>m_Range) || (flag))
		{
			m_Range=newRange;
			flag=FALSE;
		}
	}
//
	m_Scaling=1.0;
	//m_Scaling=0.008f;
	m_xTranslation = 0.0f;
	m_yTranslation = 0.0f;
	m_zTranslation = 0.0f;
	refresh();
}

void CGLKernelView::zoom(double ratio)
{
	m_Scaling*=ratio;
	if (m_Scaling<0.00001) m_Scaling=0.00001f;
	refresh();
}

void CGLKernelView::clear_tools()
{
	if (m_currentTool)
		delete m_currentTool;
	m_currentTool=NULL;
}

void CGLKernelView::set_tool(MouseTool *tool)
{
	m_currentTool=tool;
}

void CGLKernelView::screen_to_wcl(double sx, double sy, double &cx, double &cy, double &cz)
{
	GLdouble objx, objy, objz;
	double y = m_SizeY - sy;
	gluUnProject(sx, y, 0.5, modelMatrix, projMatrix, viewport, &objx, &objy, &objz); 

	cx=objx;	cy=objy;	cz=objz;
}

void CGLKernelView::wcl_to_screen(double cx, double cy, double cz, double &sx, double &sy)
{
	GLdouble winx, winy, winz;
	gluProject(cx, cy, cz, modelMatrix, projMatrix, viewport, &winx, &winy, &winz); 

	sx=winx;
	sy=m_SizeY-winy;
}

void CGLKernelView::GetUpVector(double &x, double &y, double &z)
{
	CCW_Geometry geo;
	double n[3],p1[3],p2[3],mu;
	GetViewVector(n[0],n[1],n[2]);

	screen_to_wcl(0,0,p1[0],p1[1],p1[2]);
	geo.CalPlaneLineIntersection(p1,n,n[0],n[1],n[2],0.0,mu);
	p1[0]=p1[0]+n[0]*mu;	p1[1]=p1[1]+n[1]*mu;	p1[2]=p1[2]+n[2]*mu;

	screen_to_wcl(0,-10,p2[0],p2[1],p2[2]);
	geo.CalPlaneLineIntersection(p2,n,n[0],n[1],n[2],0.0,mu);
	p2[0]=p2[0]+n[0]*mu;	p2[1]=p2[1]+n[1]*mu;	p2[2]=p2[2]+n[2]*mu;

	p2[0]=p2[0]-p1[0];	p2[1]=p2[1]-p1[1];	p2[2]=p2[2]-p1[2];
	geo.Normalize(p2);
	x=p2[0];	y=p2[1];	z=p2[2];
}

void CGLKernelView::GetViewVector(double &x, double &y, double &z)
{
	CCW_Geometry geo;
	double cx,cy,cz,d;
	double xx[3],yy[3],zz[3];

	screen_to_wcl(100,100,cx,cy,cz);
	xx[0]=cx;	yy[0]=cy;	zz[0]=cz;
	screen_to_wcl(200,200,cx,cy,cz);
	xx[1]=cx;	yy[1]=cy;	zz[1]=cz;
	screen_to_wcl(200,100,cx,cy,cz);
	xx[2]=cx;	yy[2]=cy;	zz[2]=cz;
	geo.CalPlaneEquation(x,y,z,d,xx,yy,zz);
}

void CGLKernelView::refreshRange()
{
	float newRange = 0.0;
	float oldRange = m_Range;

	for(POSITION Pos=m_displayObjList.GetHeadPosition();Pos!=NULL;)
	{
		GLEntity *entity=(GLEntity *)(m_displayObjList.GetNext(Pos));
		if (!entity) continue;
		if (!(entity->bShow)) continue;
		if ((entity->getRange())>newRange) newRange=entity->getRange();
	}

	if ((newRange>m_Range) || ((m_displayObjList.GetCount())==0))
	{
		m_Range=newRange; 
		m_Scaling=m_Scaling*(newRange/oldRange);
	}
}

void CGLKernelView::AddDisplayObj(GLEntity *entity, BOOL bRefresh)
{
	float newRange = entity->getRange();
	float oldRange = m_Range;
	
	if ((newRange>m_Range) || ((m_displayObjList.GetCount())==0))
	{
		m_Range=newRange;
		m_Scaling=m_Scaling*(newRange/oldRange); 
	}

	m_displayObjList.AddTail(entity);
	if (bRefresh) refresh();
}

UINT CGLKernelView::DisplayObjCount()
{
	return m_displayObjList.GetCount();
}

GLEntity* CGLKernelView::GetDisplayObjAt(UINT nIndex)
{
	POSITION Pos;	UINT n=1;
	for(Pos=m_displayObjList.GetHeadPosition();Pos!=NULL;)
	{
		GLEntity *tempEntity=(GLEntity *)(m_displayObjList.GetNext(Pos));
		if (n==nIndex)	return tempEntity;
		n++;
	}

	return NULL;
}

void CGLKernelView::DelDisplayObj2(GLEntity *entity)
{
	POSITION tempPos;
	if( ( tempPos = m_displayObjList.Find( entity ) ) != NULL )
		m_displayObjList.RemoveAt(tempPos);
}

void CGLKernelView::DelDisplayObj3(GLEntity *entity)
{
	float newRange,oldRange=m_Range;
	m_Range=1.0f;

	POSITION tempPos;
	if( ( tempPos = m_displayObjList.Find( entity ) ) != NULL )
		m_displayObjList.RemoveAt(tempPos);

	POSITION Pos;	BOOL flag=TRUE;
	for(Pos=m_displayObjList.GetHeadPosition();Pos!=NULL;)
	{
		GLEntity *tempEntity=(GLEntity *)(m_displayObjList.GetNext(Pos));
		newRange=tempEntity->getRange();
		if ((newRange>m_Range) || (flag))
		{
			m_Range=newRange;
			flag=FALSE;
		}
	}

	m_Scaling=m_Scaling*(m_Range/oldRange);
	refresh();	
}

void CGLKernelView::DelDisplayObj(GLEntity *entity)
{
	float newRange,oldRange=m_Range;
	m_Range=1.0f;

	POSITION tempPos;
	if( ( tempPos = m_displayObjList.Find( entity ) ) != NULL )
		m_displayObjList.RemoveAt(tempPos);

	POSITION Pos;	BOOL flag=TRUE;
	for(Pos=m_displayObjList.GetHeadPosition();Pos!=NULL;)
	{
		GLEntity *tempEntity=(GLEntity *)(m_displayObjList.GetNext(Pos));
		newRange=tempEntity->getRange();
		if ((newRange>m_Range) || (flag))
		{
			m_Range=newRange;
			flag=FALSE;
		}
	}

	m_Scaling=m_Scaling*(m_Range/oldRange);
	RedrawAll();
	//refresh();	

	delete entity;
}

void CGLKernelView::add_gllist(GLList *glList)
{
	m_glList.AddTail(glList);
}

UINT CGLKernelView::gllistCount()
{
	return m_glList.GetCount();
}

GLList* CGLKernelView::get_gllist(UINT nIndex)
{
	POSITION Pos;	UINT n=1;
	for(Pos=m_glList.GetHeadPosition();Pos!=NULL;)
	{
		GLList *glList=(GLList *)(m_glList.GetNext(Pos));
		if (n==nIndex) return glList;
		n++;
	}

	return NULL;
}

void CGLKernelView::clear_gllist()
{
	POSITION Pos;
	for(Pos=m_glList.GetHeadPosition();Pos!=NULL;)
	{
		GLList *glList=(GLList *)(m_glList.GetNext(Pos));
		if (!glList) continue;
		delete glList;
	}
	m_glList.RemoveAll();
}

void CGLKernelView::ClearDisplayObjList()
{
	m_displayObjList.RemoveAll();
	m_Range=1.0f;
	
}

void CGLKernelView::HighLightObj(GLEntity *entity)
{
	m_HighLightObj=entity;
	refresh();

	Sleep(500);

	m_HighLightObj=NULL;
	refresh();
}


void CGLKernelView::set_foreground_color(float red, float green, float blue)
{
	m_red = red;
	m_green = green;
	m_blue = blue;
}

void CGLKernelView::set_line_width(int width)
{
	m_lineWidth = width;
}

void CGLKernelView::draw_text_3d(double cx, double cy, double cz, char *s)
{
	double sx,sy;
	wcl_to_screen(cx,cy,cz,sx,sy);
	draw_text_2d((int)sx,(int)sy,s);
}

void CGLKernelView::draw_text_2d(int cx, int cy, char *s)
{
	CDC *pDC = m_pWnd->GetWindowDC();
	if (!pDC) return;

	pDC->SetTextColor(RGB((BYTE)(m_red*255),(BYTE)(m_green*255),(BYTE)(m_blue*255)));
	pDC->TextOut(cx,cy,s,(int)strlen(s));

	m_pWnd->ReleaseDC(pDC);
}

void CGLKernelView::draw_polyline_3d(UINT pointNum, position_array *pa, BOOL bFill)
{
	CDC *pDC = m_pWnd->GetWindowDC();
	if (!pDC) return;

	CPen myPen,*oldPen;
	double sx,sy;

	myPen.CreatePen(PS_SOLID,m_lineWidth,RGB((BYTE)(m_red*255),(BYTE)(m_green*255),(BYTE)(m_blue*255)));
	oldPen=pDC->SelectObject(&myPen);

	double x,y,z;

	pa->ElementAt(0,x,y,z);
	wcl_to_screen(x,y,z,sx,sy);
	pDC->MoveTo(sx,sy);
	for(UINT i=1;i<pointNum;i++)
	{
		pa->ElementAt(i,x,y,z);
		wcl_to_screen(x,y,z,sx,sy);
		pDC->LineTo((int)sx,(int)sy);
	}

	pDC->SelectObject(oldPen);

	m_pWnd->ReleaseDC(pDC);
}

void CGLKernelView::draw_polyline_3d(UINT pointNum, const float pts[], BOOL bFill)
{
	CDC *pDC = m_pWnd->GetWindowDC();
	if (!pDC) return;

	CPen myPen,*oldPen;
	double sx,sy;

	myPen.CreatePen(PS_SOLID,m_lineWidth,RGB((BYTE)(m_red*255),(BYTE)(m_green*255),(BYTE)(m_blue*255)));
	oldPen=pDC->SelectObject(&myPen);
	
	wcl_to_screen(pts[0],pts[1],pts[2],sx,sy);
	pDC->MoveTo((int)sx,(int)sy);
	for(UINT i=1;i<pointNum;i++)
	{
		wcl_to_screen(pts[i*3+0],pts[i*3+1],pts[i*3+2],sx,sy);
		pDC->LineTo((int)sx,(int)sy);
	}

	pDC->SelectObject(oldPen);

	m_pWnd->ReleaseDC(pDC);
}

void CGLKernelView::draw_polyline_2d(UINT pointNum, const float pts[], BOOL bFill)
{
	CDC *pDC = m_pWnd->GetWindowDC();
	if (!pDC) return;

	CPen myPen,*oldPen;
	int sx,sy;

	pDC->SetROP2(R2_XORPEN);
	myPen.CreatePen(PS_SOLID,m_lineWidth,RGB((BYTE)(m_red*255),(BYTE)(m_green*255),(BYTE)(m_blue*255)));
	oldPen=pDC->SelectObject(&myPen);
	sx=(int)pts[0];	sy=(int)pts[1];
	pDC->MoveTo(sx,sy);
	for(UINT i=1;i<pointNum;i++)
	{
		sx=(int)pts[i*2+0];	sy=(int)pts[i*2+1];
		pDC->LineTo(sx,sy);
	}
	pDC->SelectObject(oldPen);
	pDC->SetROP2(R2_COPYPEN);

	m_pWnd->ReleaseDC(pDC);
}

BOOL CGLKernelView::pick_entity(const pick_event& pe, UINT entityType, GLEntity*& ent)
{
	GLuint selectBuf[BUFSIZE];
	GLint viewport[4];
	GLint hits;
	int x,y;

	CDC *pDC = m_pWnd->GetWindowDC();
	if (!pDC) return FALSE;
	wglMakeCurrent(pDC->m_hDC,m_hRC);

	x=pe.x;	y=pe.y;
	glSelectBuffer(BUFSIZE,selectBuf);
	glGetIntegerv(GL_VIEWPORT,viewport);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();

	glRenderMode(GL_SELECT);
	glLoadIdentity();
	gluPickMatrix((GLdouble)x,(GLdouble)(viewport[3]-y),10.0,10.0,viewport);

	GLPickInit();
	GLPickDraw(entityType);

	hits=glRenderMode(GL_RENDER);

	//////////////////////////////////////////////////////////////////////////
	//	The following lines analyze the hit information
	GLuint *ptr=(GLuint*)selectBuf;
	GLuint names;
	UINT currentName=0,z1,z2,zmin,selectName,selectZ;

	for(int i=0;i<hits;i++)
	{
		names=*ptr;		ptr++;
		z1=*ptr;		ptr++;
		z2=*ptr;		ptr++;
		for(UINT j=0;j<names;j++)
		{
			currentName=*ptr;	ptr++;
		}
		if (z1<z2)
			zmin=z1;
		else
			zmin=z2;

		if (i==0)
		{
			selectName=currentName;
			selectZ=zmin;
		}
		else
		{
			if (zmin<selectZ)
			{
				selectName=currentName;
				selectZ=zmin;
			}
		}
	}

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
    wglMakeCurrent(NULL,NULL);

	m_pWnd->ReleaseDC(pDC);

	if (currentName==0) return FALSE;

	ent=GetPickObjByName(entityType,selectName);

	if (!ent) return FALSE;

	return TRUE;
}

void CGLKernelView::renderBitmapCharacher(float x, float y, float z, void *font,char *string)
{
	font=GLUT_BITMAP_HELVETICA_18;
  char *c;
  glRasterPos3f(x, y,z);
  for (c=string; *c != '\0'; c++) {
    glutBitmapCharacter(font, *c);
  }
}