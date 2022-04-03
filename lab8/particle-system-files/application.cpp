#include "application.h"

#include <iostream>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>

#ifndef __APPLE__
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#endif

using namespace std;
enum { NONE, AMBIENT, DIFFUSE, SPECULAR, NUM_MODES };

// Particle struct 
struct Particle
{
    float m = 0.0; // for the mass of a partcle
    vec3 p; // for the positive of a particle
    // below are varibale have direction,so its vector
    vec3 v; // for the velocity of a particle
    vec3 f; // for the force applied on a particle
    vec3 c; // for the color of the particle
    float d; // for duration of the particle

    // memeber fucntion
    // update v and x(p) with an Forward Euler Step
    void Euler_Step(float h)
    {
        // position step (update position)
        p += h * v;
        // velocity step (update velocity)
        v += (h / m) * f;
    }

    // reset force to 0 vector;
    void Reset_Forces() {
        f = vec3(0, 0, 0); 
    }
    // for handle the partice hit the ground
    void Handle_Collision(float damping, float coeff_resititution) {
        // for restitution
        // if x_y < 0 (x_x,x_y,x_z) = (p[0], p[1], p[2])
        if(p[1] < 0)
        {
            //reset x_y to 0
            p[1] = 0;
            // in additon to restitution
            // we apply damping forces
            // if v_y < 0
            if(v[1] < 0)
            {
                // updated our velocity v_y
                // bounce upward after hit the ground
                v[1] = -coeff_resititution * v[1];
                //upated our v_x, and v_z
                v[0] = damping * v[0];
                v[2] = damping * v[2];
            }
        }
    }
};
// give random float depend on range low nad high
float random(float low, float high) {
    // it will generate a number from range low to high number
    float r = low + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(high-low)));
    // returns a random float between low and high
    return r;
}

// global variable(s) to keep a list of particles
vector<Particle> list_p;
// helper functions for add randomly initialized particles
// this fucntion generates n random particles, and appends to the particle vector.
// add new partcles and initilze
void Add_Particles(int n) {
    for(int i = 0; i < n; i++) {
        //inclusive generate a number from 0.0 to 1.0
        Particle part;
        // particle inital values that suggest by lab 
        // mass of particle: 1
        part.m = 1.0f;
        // start position of a particle, x: (random(-0.2,0.2), 0.05,random(-0.2,0.2))
        part.p = vec3(random(-0.2f,0.2f), 0.05, random(-0.2f,0.2f));
        // v: (10*x.x,random(1,10),10*x.z)
        part.v = vec3(10 * part.p[0], random(1.0f, 10.0f), 10 * part.p[2]);
        // color of the particle: yellow = (255,255,0) or (1,1,0)
        part.c = vec3(1,1,0);
        // duration inital value at 0
        part.d = 0.0f;
        list_p.push_back(part); // store particles to our list
    }
}

//global helper function that returns the interpolated color
vec3 Get_Particle_Color(float d) {
    //if d < 0.1: return yellow
    if (d < 0.1) {
        return vec3(1,1,0);
    }
    //else if d < 1.5: return an interpolated value from yellow to red.
    else if (d < 1.5) {
        // interpolated value from yellow to red
        // yellow is (1,1,0);
        // red is (1,0,0)
        // interpolate result is (color2 - color1) * fraction + color1
        float fraction = (d - 0.1)/(1.5-0.1); 
        float R = (1-1) * fraction + 1;
        float G = (0-1) * fraction + 1;
        float B = (0-0) * fraction + 0;
        return vec3(R,G,B);
    }
    //else if d < 2: return red.
    else if (d < 2.0) {
        return vec3(1,0,0);
    }
    //else if d < 3: return an interpolated value from red to grey.
    else if (d < 3.0) {
        // red (1,0,0)
        // grey (0.5,0.5,0.5)
        // interpolate result is (color2 - color1) * fraction + color1
        float fraction = (d - 2.0)/(3.0-2.0);
        float R = (0.5-1) * fraction + 1;
        float G = (0.5-0) * fraction + 0;
        float B = (0.5-0) * fraction + 0;
        return vec3(R,G,B);
    }
    //else return grey.
    else {
        return vec3(0.5,0.5,0.5);
    }
} 

void draw_grid(int dim);
void draw_obj(obj *o, const gl_image_texture_map& textures);

void set_pixel(int x, int y, float col[3])
{
    // write a 1x1 block of pixels of color col to framebuffer
    // coordinates (x, y)
    //glRasterPos2i(x, y);
    //glDrawPixels(1, 1, GL_RGB, GL_FLOAT, col);

    // use glVertex instead of glDrawPixels (faster)
    glBegin(GL_POINTS);
    glColor3fv(col);
    glVertex2f(x, y);
    glEnd();
}

application::application()
    : raytrace(false), rendmode(SPECULAR), paused(false), sim_t(0.0),draw_volcano(true),h(0.015)
{
}

application::~application()
{
}

// triggered once after the OpenGL context is initialized
void application::init_event()
{

    cout << "CAMERA CONTROLS: \n  LMB: Rotate \n  MMB: Move \n  RMB: Zoom" << endl;
    cout << "KEYBOARD CONTROLS: \n";
    cout << "  'p': Pause simulation\n";
    cout << "  'v': Toggle draw volcano" << endl;

    const GLfloat ambient[] = { 0.0, 0.0, 0.0, 1.0 };
    const GLfloat diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
    const GLfloat specular[] = { 1.0, 1.0, 1.0, 1.0 };

    // enable a light
    glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT1, GL_SPECULAR, specular);
    glEnable(GL_LIGHT1);

    // set global ambient lighting
    GLfloat amb[] = { 0.4, 0.4, 0.4, 1.0 };
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, amb);

    // enable depth-testing, colored materials, and lighting
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glDisable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHTING);

    // normalize normals so lighting calculations are correct
    // when using GLUT primitives
    glEnable(GL_RESCALE_NORMAL);

    // enable smooth shading
    glShadeModel(GL_SMOOTH);

    glClearColor(0,0,0,0);

    set_camera_for_box(vec3(-3,-2,-3),vec3(3,5,3));

    t.reset();
    o.load("crater.obj");

    // loads up all the textures referenced by the .mtl file
    const std::map<std::string, obj::material>& mats = o.get_materials();
    for (std::map<std::string, obj::material>::const_iterator i = mats.begin();
        i != mats.end(); ++i
        )
    {
        if (!i->second.map_kd.empty()) {
            string filename = i->second.map_kd;

            // add texture if we have not already loaded it
            if (texs.find(filename) == texs.end()) {
                gl_image_texture *tex = new gl_image_texture();
                if (tex->load_texture(filename) != SUCCESS) {
                    cout << "could not load texture file: " << filename << endl;
                    exit(0);
                }
                texs[filename] = tex;
            }
        }
    }
    //Create 10 new particles in the init event
    Add_Particles(10);
}

// triggered each time the application needs to redraw
void application::draw_event()
{
    apply_gl_transform();

    const GLfloat light_pos1[] = { 0.0, 10.0, 0.0, 1 };
    glLightfv(GL_LIGHT1, GL_POSITION, light_pos1);

    if (!paused) {
        //
        //ADD NEW PARTICLES
        //
        //
        // SIMULATE YOUR PARTICLE HERE.
        //
        //
        //
        // UPDATE THE COLOR OF THE PARTICLE DYNAMICALLY
        //
        // add new particles
        Add_Particles(20);
        // simluate our parteicle
        for(auto& this_part : list_p)
        {
          
          // updated our duration value before we draw the particles
          // updated with time-step h: d = d + h 
          this_part.d = this_part.d + h;

          // update velocity(v) and position (x or p)
          this_part.Euler_Step(h);
          // reset force
          this_part.Reset_Forces();
          // add force
          this_part.f = this_part.f  + vec3(0, -1.0 * this_part.m * 9.8, 0);
          // handle the collision when particle hit the ground
          this_part.Handle_Collision(.5, .5);
          //this_part.Handle_Collision(1, 1);
          //this_part.Handle_Collision(0, 0);
          //this_part.Handle_Collision(0.5, 1);

          // update the particle's color
          this_part.c = Get_Particle_Color(this_part.d);
        }

    }

    glLineWidth(2.0);
    glEnable(GL_COLOR_MATERIAL);
    glBegin(GL_LINES);
        //
        //
        // DRAW YOUR PARTICLE USING GL_LINES HERE
        //
        // glVertex3f(...) endpoint 1
        // glVertex3f(...) endpoint 2
        //
        //
        //
        // Draw each particle p as a line from p.x to p.x+0.04*p.v, and with color
        for(auto& this_part : list_p) {
            // color
            glColor3f(this_part.c[0],this_part.c[1],this_part.c[2]);
            // position
            glVertex3f(
                this_part.p[0],    // x
                this_part.p[1],    // y
                this_part.p[2]);   // z
            // velocity
            glVertex3f(
                this_part.p[0] + (0.04 * this_part.v[0]),  // x + s*v_x
                this_part.p[1] + (0.04 * this_part.v[1]),  // y + s*v_y
                this_part.p[2] + (0.04 * this_part.v[2])); // z + s*v_z
        }
    glEnd();

    // draw the volcano
    if(draw_volcano){
        glEnable(GL_LIGHTING);
        glPushMatrix();
        glScalef(0.2,0.3,0.2);
        draw_obj(&o, texs);
        glPopMatrix();
        glDisable(GL_LIGHTING);
    }


    glColor3f(0.15, 0.15, 0.15);
    draw_grid(40);

    //
    // This makes sure that the frame rate is locked to close to 60 fps.
    // For each call to draw_event you will want to run your integrate for 0.015s
    // worth of time.
    //
    float elap = t.elapsed();
    if (elap < h) {
        usleep(1e6*(h-elap));
		//sleep(h-elap);
    }
    t.reset();
}

// triggered when mouse is clicked
void application::mouse_click_event(int button, int button_state, int x, int y)
{
}

// triggered when mouse button is held down and the mouse is
// moved
void application::mouse_move_event(int x, int y)
{
}

// triggered when a key is pressed on the keyboard
void application::keyboard_event(unsigned char key, int x, int y)
{

    if (key == 'r') {
        sim_t = 0;
    } else if (key == 'p') {
        paused = !paused;
    } else if (key == 'v') {
        draw_volcano=!draw_volcano;
    } else if (key == 'q') {
        exit(0);
    }
}

void draw_grid(int dim)
{
    glLineWidth(2.0);


    //
    // Draws a grid along the x-z plane
    //
    glLineWidth(1.0);
    glBegin(GL_LINES);

    int ncells = dim;
    int ncells2 = ncells/2;

    for (int i= 0; i <= ncells; i++)
    {
        int k = -ncells2;
        k +=i;
        glVertex3f(ncells2,0,k);
        glVertex3f(-ncells2,0,k);
        glVertex3f(k,0,ncells2);
        glVertex3f(k,0,-ncells2);
    }
    glEnd();

    //
    // Draws the coordinate frame at origin
    //
    glPushMatrix();
    glScalef(1.0, 1.0, 1.0);
    glBegin(GL_LINES);

    // x-axis
    glColor3f(1.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(1.0, 0.0, 0.0);

    // y-axis
    glColor3f(0.0, 1.0, 0.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 1.0, 0.0);

    // z-axis
    glColor3f(0.0, 0.0, 1.0);
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, 1.0);
    glEnd();
    glPopMatrix();
}

void draw_obj(obj *o, const gl_image_texture_map& textures)
{
    glDisable(GL_COLOR_MATERIAL);

    // draw each polygon of the mesh
    size_t nfaces = o->get_face_count();
    for (size_t i = 0; i < nfaces; ++i)
    {
        const obj::face& f = o->get_face(i);

        // sets the material properties of the face
        const obj::material& mat = o->get_material(f.mat);
        if (!mat.map_kd.empty()) {
            gl_image_texture_map::const_iterator it = textures.find(mat.map_kd);
            if (it != textures.end()) {
                gl_image_texture* tex = it->second;
                tex->bind();
            }
            GLfloat mat_amb[] = { 1, 1, 1, 1 };
            GLfloat mat_dif[] = { mat.kd[0], mat.kd[1], mat.kd[2], 1 };
            GLfloat mat_spec[] = { mat.ks[0], mat.ks[1], mat.ks[2], 1 };
            glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_dif);
            glMaterialfv(GL_FRONT, GL_SPECULAR, mat_spec);
        } else {
            GLfloat mat_amb[] = { mat.ka[0], mat.ka[1], mat.ka[2], 1 };
            GLfloat mat_dif[] = { mat.kd[0], mat.kd[1], mat.kd[2], 1 };
            GLfloat mat_spec[] = { mat.ks[0], mat.ks[1], mat.ks[2], 1 };
            glMaterialfv(GL_FRONT, GL_AMBIENT, mat_amb);
            glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_dif);
            glMaterialfv(GL_FRONT, GL_SPECULAR, mat_spec);
        }
        glMaterialf(GL_FRONT, GL_SHININESS, mat.ns);

        if (!glIsEnabled(GL_TEXTURE_2D)) glEnable(GL_TEXTURE_2D);

        // draws a single polygon
        glBegin(GL_POLYGON);
        for (size_t j = 0; j < f.vind.size(); ++j)
        {
            // vertex normal
            if (f.nind.size() == f.vind.size()) {
                const float *norm = o->get_normal(f.nind[j]);
                glNormal3fv(norm);
            }

            // vertex UV coordinate
            if (f.tex.size() > 0) {
                const float* tex = o->get_texture_indices(f.tex[j]);
                glTexCoord2fv(tex);
            }

            // vertex coordinates
            const float *vert = o->get_vertex(f.vind[j]);
            glVertex3fv(vert);
        }
        glEnd();
    }
    glDisable(GL_TEXTURE_2D);
    glEnable(GL_COLOR_MATERIAL);
}
