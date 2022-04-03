#include "sphere.h"
#include "ray.h"

// Determine if the ray intersects with the sphere
Hit Sphere::Intersection(const Ray& ray, int part) const
{
    //TODO;
    Hit result;
    result.object = NULL;
    // ray-sphere-intersection: (d · d)t^2 + 2d · (e − c)t + (e − c) · (e − c) − R^2 = 0.(at^2+bt+c=0)
    // (e-c) is endpoint- center
    vec3 ray_to_sphere = ray.endpoint - this->center; 
    float a = dot(ray.direction, ray.direction);  // (d · d)
    float b = 2 * dot(ray.direction, ray_to_sphere); // 2d · (e − c)
    float c = dot(ray_to_sphere, ray_to_sphere) - (radius * radius); // (e − c)t + (e − c) · (e − c) − R^2

    // quadratic formula : [-b ± √(b² - 4ac)]/2a
    // our discriminant determine by (b² - 4ac)
    float discriminant = (b*b) - 4 * a * c;

    // if discrimiinant large than 0 indicate we have solutions
    // Other like '< 0' is no solution, and '= 1' is 1 solution (ray tracer is treat as no solution)
    if (discriminant > 0) { // if discriminant large than zero (positive)
        // we have two solution t1 and t2
        // t1 = [-b + √(b² - 4ac)]/2a
        float t_1 = (-b + sqrt(discriminant))/(2 * a);
        // t2 = [-b - √(b² - 4ac)]/2a
        float t_2 = (-b - sqrt(discriminant))/(2 * a);
        if (t_1 >= small_t && t_2 >= small_t) {  // if both t_1 and t_2 are larage than small_t
            // calculate the smallest value, and store it
            float t = std::min(t_1, t_2);
			result.object = this;
			result.dist = t;
		}
        else if (t_1 >= small_t && t_2 < small_t) { // if t_1 >= small_t and t_2 < small_T (if t_1 > 0)
            result.object = this;
			result.dist = t_1;
        }
        else if (t_1 < small_t && t_2 >= small_t) { // if t_1 < small_t and t_2 >= small_T (if t_2 > 0)
            result.object = this;
			result.dist = t_2;
        }
		else {
			result.object = NULL;
		}
    }
    return result;
}

vec3 Sphere::Normal(const vec3& point, int part) const
{
    vec3 normal;
    //TODO; // compute the normal direction
    // (endpoint - center).normalized 
    normal = (point - center).normalized();
    return normal;
}

Box Sphere::Bounding_Box(int part) const
{
    Box box;
    TODO; // calculate bounding box
    return box;
}
