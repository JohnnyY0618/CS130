#include "render_world.h"
#include "flat_shader.h"
#include "object.h"
#include "light.h"
#include "ray.h"

extern bool disable_hierarchy;

Render_World::Render_World()
    :background_shader(0),ambient_intensity(0),enable_shadows(true),
    recursion_depth_limit(3)
{}

Render_World::~Render_World()
{
    delete background_shader;
    for(size_t i=0;i<objects.size();i++) delete objects[i];
    for(size_t i=0;i<lights.size();i++) delete lights[i];
}

// Find and return the Hit structure for the closest intersection.  Be careful
// to ensure that hit.dist>=small_t.
Hit Render_World::Closest_Intersection(const Ray& ray)
{
    //TODO;
        // use to store ray object
    Hit closest_hit;
    // since we have not find the ray object intersection
    // initial it to NULL
    closest_hit.object = NULL;
    float min_t =std::numeric_limits<float>::max();
    for (unsigned i = 0; i < objects.size(); i++) {
        // use temp variable(Hit) to store the closest hit with the object
        Hit temp_hit = objects[i]->Intersection(ray, -1);
        // if Hit is the closest so far and larger than small t
        if (temp_hit.dist < min_t && temp_hit.dist > small_t) {
            // Store the hit as the closest hit
            closest_hit = temp_hit;
            min_t = temp_hit.dist;
        }
    }
    return closest_hit;
    //return {};
}

// set up the initial view ray and call
void Render_World::Render_Pixel(const ivec2& pixel_index)
{
    //TODO; // set up the initial view ray here
    Ray ray;
    // s − e is the ray’s direction.  s is point of the image plane,  e (end point) is the camera position (from camera class).
    //direction is a unit vector from the camera position to the world position of the pixel. (so use normalized function)
    vec3 direction = (camera.World_Position(pixel_index) - camera.position).normalized();
    ray.endpoint =  camera.position;
    ray.direction = direction;
    
    vec3 color=Cast_Ray(ray,1);
    camera.Set_Pixel(pixel_index,Pixel_Color(color));
}

void Render_World::Render()
{
    if(!disable_hierarchy)
        Initialize_Hierarchy();

    for(int j=0;j<camera.number_pixels[1];j++)
        for(int i=0;i<camera.number_pixels[0];i++)
            Render_Pixel(ivec2(i,j));
}

// cast ray and return the color of the closest intersected surface point,
// or the background color if there is no object intersection
vec3 Render_World::Cast_Ray(const Ray& ray,int recursion_depth)
{
    vec3 color;
    //TODO; // determine the color here
    // Get the closest hit with an object
    Hit hit = Closest_Intersection(ray);
    //  If there is an intersection set color using the object Shade Surface
    //          Shade Surface receives as parameters: ray, intersection point, normal at 
    //          the intersection point and recursion depth. 
    //          object is null represent no intersection, otherwise, it has intersection.
    if (hit.object != NULL) { // if there is an intersection
        color = hit.object->material_shader->Shade_Surface(ray, ray.Point(hit.dist), hit.object->Normal(ray.Point(hit.dist),hit.part), recursion_depth);

    }
    else { // if there is no intersection
        color = background_shader->Shade_Surface(ray, vec3(0,0,0), vec3(0,0,0), recursion_depth);
    }
    return color;
}

void Render_World::Initialize_Hierarchy()
{
    TODO; // Fill in hierarchy.entries; there should be one entry for
    // each part of each object.

    hierarchy.Reorder_Entries();
    hierarchy.Build_Tree();
}
