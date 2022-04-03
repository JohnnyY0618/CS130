#include "plane.h"
#include "ray.h"
#include <cfloat>
#include <limits>

// Intersect with the half space defined by the plane.  The plane's normal
// points outside.  If the ray starts on the "inside" side of the plane, be sure
// to record a hit with t=0 as the first entry in hits.
Hit Plane::Intersection(const Ray& ray, int part) const
{
    //TODO;
    float t;
    float denominator;
    Hit intersect = {NULL,0,0};
    // t = ((p-e)*n)/(u*n)
    t = dot(x1-ray.endpoint, normal)/dot(ray.direction,normal);
    // check if our denominator is zero or not
    denominator = dot(ray.direction,normal);
    if (!denominator) { // if denominator is zero
        return intersect; // return NULL
    }
    else {
        if (t < small_t) {  // check our t boundary, if  t < small_t
            return intersect; // return NULL
        }
        else {
            intersect.object = this;
            intersect.dist = t;
        }
        

    }
    return intersect;
}

vec3 Plane::Normal(const vec3& point, int part) const
{
    return normal;
}

// There is not a good answer for the bounding box of an infinite object.
// The safe thing to do is to return a box that contains everything.
Box Plane::Bounding_Box(int part) const
{
    Box b;
    b.hi.fill(std::numeric_limits<double>::max());
    b.lo=-b.hi;
    return b;
}
