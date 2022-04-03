#include "mesh.h"
#include <fstream>
#include <string>
#include <limits>

#include "plane.h"

// Consider a triangle to intersect a ray if the ray intersects the plane of the
// triangle with barycentric weights in [-weight_tolerance, 1+weight_tolerance]
static const double weight_tolerance = 1e-4;

// Read in a mesh from an obj file.  Populates the bounding box and registers
// one part per triangle (by setting number_parts).
void Mesh::Read_Obj(const char* file)
{
    std::ifstream fin(file);
    if(!fin)
    {
        exit(EXIT_FAILURE);
    }
    std::string line;
    ivec3 e;
    vec3 v;
    box.Make_Empty();
    while(fin)
    {
        getline(fin,line);

        if(sscanf(line.c_str(), "v %lg %lg %lg", &v[0], &v[1], &v[2]) == 3)
        {
            vertices.push_back(v);
            box.Include_Point(v);
        }

        if(sscanf(line.c_str(), "f %d %d %d", &e[0], &e[1], &e[2]) == 3)
        {
            for(int i=0;i<3;i++) e[i]--;
            triangles.push_back(e);
        }
    }
    number_parts=triangles.size();
}

// Check for an intersection against the ray.  See the base class for details.
Hit Mesh::Intersection(const Ray& ray, int part) const
{
    //TODO;
    Hit result;
    result.object = NULL;
    if (part >= 0) { 
        if (Intersect_Triangle(ray, part, result.dist)) {
            result.object = this;
            result.part = part;
        }
    }
    else {
        result.dist = std::numeric_limits<float>::max();// assume our result to max distance
        double temp_dist; // for store temppray distance
        for (unsigned i = 0; i < triangles.size(); i++) {
            if (Intersect_Triangle(ray, i, temp_dist)) {
                if (temp_dist < result.dist) {
                    result.object = this;
                    result.dist = temp_dist;
                    result.part = i;
                }
            }
        }
    }
    return result;
}

// Compute the normal direction for the triangle with index part.
vec3 Mesh::Normal(const vec3& point, int part) const
{
    assert(part>=0);
    //TODO;
    // get a current trinagle vertex repersent a, b ,c 
    ivec3 current_triangle  = triangles[part];
    vec3 vertex_a = vertices[current_triangle[0]];
	vec3 vertex_b = vertices[current_triangle[1]];
	vec3 vertex_c = vertices[current_triangle[2]];
    // The normal = u x v
    // u = the vector that between vertex b and vertex a
    vec3 u = vertex_b - vertex_a;
    // v = the vector that bewteen vertex c and vertex a 
    vec3 v = vertex_c - vertex_a;
    vec3 normal = cross(u, v).normalized();

    return normal;
    //return vec3();
}

// This is a helper routine whose purpose is to simplify the implementation
// of the Intersection routine.  It should test for an intersection between
// the ray and the triangle with index tri.  If an intersection exists,
// record the distance and return true.  Otherwise, return false.
// This intersection should be computed by determining the intersection of
// the ray and the plane of the triangle.  From this, determine (1) where
// along the ray the intersection point occurs (dist) and (2) the barycentric
// coordinates within the triangle where the intersection occurs.  The
// triangle intersects the ray if dist>small_t and the barycentric weights are
// larger than -weight_tolerance.  The use of small_t avoid the self-shadowing
// bug, and the use of weight_tolerance prevents rays from passing in between
// two triangles.
bool Mesh::Intersect_Triangle(const Ray& ray, int tri, double& dist) const
{
    //TODO;
    vec3 vertex_a = vertices[triangles[tri][0]];
    vec3 vertex_b = vertices[triangles[tri][1]];
    vec3 vertex_c = vertices[triangles[tri][2]];

    // check if object hit or not
    Hit hit = Plane(vertex_a, Normal(vertex_a, tri)).Intersection(ray, tri);
    if (!hit.object) {
        return false;
    }

    // vectors for helping calculation
    vec3 u = ray.direction;
    vec3 v = vertex_b - vertex_a;
    vec3 w = vertex_c - vertex_a;
    vec3 y = ray.endpoint - vertex_a;
    // deonminator for the weight
    float deonminator = dot(cross(u,v), w);
    if (!deonminator) { // divide zero is no allow and object is not hit
        return false;
    }
    else {
        // barycentric weights calculation
        float beta = dot(cross(w, u), y) / deonminator;
        float gamma = dot(cross(u, v), y) / deonminator;
        float alpha = 1 - beta - gamma;

        if (alpha > weight_tolerance && beta > weight_tolerance && gamma > weight_tolerance) {
            dist = hit.dist;
            return true;
        }
    }
    return false;
}

// Compute the bounding box.  Return the bounding box of only the triangle whose
// index is part.
Box Mesh::Bounding_Box(int part) const
{
    Box b;
    TODO;
    return b;
}
