// Name: Jianeng Yang
// Quarter, Year: winter, 2022
// Lab: 7
//
// This file is to be modified by the student.
// main.cpp
////////////////////////////////////////////////////////////
#ifndef __APPLE__
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#endif

#include <vector>
#include <cstdio>
#include <math.h>
#include "vec.h"
#include <iostream>

using namespace std;
const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;

// global vector to store the control points
// vector 2 for x and y
std::vector<vec2> control_points;

float factorial(int n)
{
    // 0! = 1, 1! = 1, so return 1
    if (n <= 1) {
        return 1;
    }
    // recursive call
    // then the final will return 
    // n! = n*(n-1)*(n-2)*...3*2*1
    return n * factorial(n-1);
}

float combination(int n, int i)
{
    // combination formula of nCi is n!/((n – i)!*i!) 
    // The numerator and denominator are all factorial, then we can use the factorial function we create
    float c = 0;
    c = factorial(n)/( factorial(n-i)*factorial(i) );
    return c;
  
}

float binomial(int n, int i, float t)
{
    //float bc = 0.0; // for store the bezier curve
    //int c= 0; // for store the combination value
    //c = combination(n,i);
    //for (int index = 0; index < n; index++) {
        // Bezier curve = nCi * t^i *(1-t)^(n-i)* P_i
    //    bc += c*pow(t, index)*pow((1-t), n-index) * points[index];
   // }
    //return bc;
    return combination(n, i) * pow(t, i) * pow((1 - t), (n - i));
  
}

void GL_render()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glutSwapBuffers();

    //glBegin(GL_LINES);
    glBegin(GL_LINE_STRIP);
    glColor3f(1.0f,0.0f,0.0f);

    // TODO: just for example, you will need to change this.
    //glVertex2f(-.5f,-.5f);
    //glVertex2f(.5f,-.5f);
    //glVertex2f(.5f,.5f);
    //glVertex2f(-.5f,.5f);
    // to drow line segments we need at least two points
    if (control_points.size() >= 2) {
        // iterate t between 0 and 1 in steps of 0.01
        for (float t = 0.01; t < 1.0; t+= 0.01) {
            vec2 bc;
            for (size_t i = 0; i < control_points.size(); i++) {
                // Bezier curve = nCi * t^i *(1-t)^(n-i)* P_i
                bc += binomial(control_points.size(), i, t) * control_points[i];
            }
            glVertex2f(bc[0], bc[1]);
        }
    }

    glEnd();
    glFlush();
}

void GL_mouse(int button,int state,int x,int y)
{
    y=WINDOW_HEIGHT-y;
    GLdouble mv_mat[16];
    GLdouble proj_mat[16];
    GLint vp_mat[4];
    glGetDoublev(GL_MODELVIEW_MATRIX,mv_mat);
    glGetDoublev(GL_PROJECTION_MATRIX,proj_mat);
    glGetIntegerv(GL_VIEWPORT,vp_mat);

    if(button==GLUT_LEFT_BUTTON && state==GLUT_DOWN){
        double px,py,dummy_z; // we don't care about the z-value but need something to pass in
        gluUnProject(x,y,0,mv_mat,proj_mat,vp_mat,&px,&py,&dummy_z);
        // TODO: the mouse click coordinates are (px,py).
        control_points.push_back(vec2(px,py));
        glutPostRedisplay();
    }
}

//Initializes OpenGL attributes
void GLInit(int* argc, char** argv)
{
    glutInit(argc, argv);
    //glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);

    //glMatrixMode(GL_PROJECTION_MATRIX);
    //glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, -1, 1);
    glutCreateWindow("CS 130 - <Insert Name Here>");
    glutDisplayFunc(GL_render);
    glutMouseFunc(GL_mouse);
}

int main(int argc, char** argv)
{
    GLInit(&argc, argv);
    glutMainLoop();
    return 0;
}
