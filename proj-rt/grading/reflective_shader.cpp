#include "reflective_shader.h"
#include "ray.h"
#include "render_world.h"

vec3 Reflective_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& normal,int recursion_depth) const
{
    vec3 color;
    //TODO; // determine the color
    // call shade_surface funciton in phong_sader to get our color(phong reflection)
    color = shader->Shade_Surface(ray, intersection_point, normal, recursion_depth);
    // if our recursion depth large than recursion depth limit
    if (recursion_depth >= world.recursion_depth_limit) {
        color *= (1 - reflectivity);  // set color to 1 - reflectivity
    }
    else { // else recursion depth is small than recusion depth limit
        vec3 d = ray.direction; // get ray dirction
        // r = d − 2(d · n)n; 
        Ray reflect_ray(intersection_point, d - 2 * dot(d, normal) * normal); // reflection ray
        vec3 reflect_color = world.Cast_Ray(reflect_ray, recursion_depth+1);
        color = (1 - reflectivity) * color + reflectivity * reflect_color;
    }
    return color;
}
