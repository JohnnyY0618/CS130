#include "light.h"
#include "phong_shader.h"
#include "ray.h"
#include "render_world.h"
#include "object.h"

vec3 Phong_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& normal,int recursion_depth) const
{
    vec3 color; // return variable combine with ambinent, diffuse, specular
    //TODO; //determine the color
    vec3 ambient; // for store the ambinet of the world/scense
    vec3 diffuse; // for store diffuse of all the lights
    vec3 specular;  // for store specfular of the lights
    vec3 l; //for store light sources

    float diffuse_intensity; // for store diffuse intensity
    float specular_intensity; // for story specular intensity

    vec3 r; // for store reflection direction use to cacluate specular
    

    // ambient(I_a) = R_a * L_a
    // R_a is color of the object
    // L_a is the intensit of the ambinet light around the scense/world
    // ambient is conbination of: world.ambient_color, world.ambienta_intensity, and color_ambient
    ambient = color_ambient * world.ambient_color * world.ambient_intensity;

    // use loop to calculate diffuse, specular for all lights(light source) in the scense
    for (unsigned int i = 0; i < world.lights.size(); i++) {
        // get the light vector that touch the scene
        l = world.lights[i]->position - intersection_point;
        
        
        // check the object is in shaw or not
        // if world.enable_shows is true
        if (world.enable_shadows) { 
            // then check if there is an object between intersectin point and the light source
            // ray object: intersectin point to the light soucrce
            Ray intersect_point_to_light(intersection_point, l);
            Hit hit = world.Closest_Intersection(intersect_point_to_light);
            // if there is hit object
            if (hit.object != NULL) {
                // if the object is blocking all the light sources
                if ((intersect_point_to_light.Point(hit.dist)-intersection_point).magnitude() < l.magnitude()) {
                    continue;  // got to next light source, which we only use ambient for current light source and skip to calculate the current light source's diffuse and specular
                }
            }
        }
        // (Diffuse)I_d = R_d*L_d*max(n*l,0)
        // R_d is the color of the object that is reflected off
        // (max(n*l,0)) is intensity of the light that hiting the object in the scene/world
        diffuse_intensity =  (std::max(dot(normal.normalized(), l.normalized()), 0.0));
        // note: use Emittedd_ligth fucntio nto get ligth color at proper intentsity
        // find and add each light's diffuse
        diffuse += color_diffuse * world.lights[i]->Emitted_Light(l) * diffuse_intensity;

        // (Specular)I_s = R_s*L_s*max(cos phi,0)^(alpha)
        // R_s is color of the object that reflected of as shiny spot
        // L_s is intensity of the light that hiting the object in the scene/world
        // calculate the light reflected direction: r = (2*(l dot n)*n-l) (l and r are normalized)
        r = 2 * ( dot(l.normalized(), normal.normalized()) ) * normal.normalized() - l.normalized();
        // specluar intensity is max(r dot c, 0)^(specular_power) or (max(cos phi,0)^(alpha));
        // alpha is specular power constant
        specular_intensity = std::pow( std::max(dot(r, -(ray.direction)), 0.0), specular_power);
        // find and add each light's specular 
        specular += color_specular * world.lights[i]->Emitted_Light(l) * specular_intensity;
    }
    // phong reflection
    color = ambient + diffuse + specular;
    return color;
}
