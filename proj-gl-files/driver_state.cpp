#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=0;
    state.image_depth=0;
    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;

    // initilze our color and depth of the image 
    // image size is width*height
    state.image_color = new pixel[width * height];
    state.image_depth = new float[width * height];
    // loop for initlize every pixel in our image
    for(int i = 0; i < width * height; i++)
    {
        state.image_color[i] = make_pixel(0, 0, 0);  // initlize to black
        state.image_depth[i] = MAX_FLOATS_PER_VERTEX; // max float per vertex can hold
    }
}

// This function will be called to render the data that has been stored in this class.
// Valid signs of type are:
//   render_type::triangle - Eac_edgeh group of three vertices corresponds to a triangle.
//   render_type::indexed -  Eac_edgeh group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    //std::cout<<"TODO: implement rendering."<<std::endl;
    // store render type of tiangle 
    data_geometry tri_output[3];  // for store the geo vector that we send to clip triangle
    data_geometry temporary_geo[3];  // for store temporary geo
    data_vertex input_geo[3];  // fore store our input geo for vertex_shader
    int num_triangles = 0; // for store num of triangles
    int vertex_index = 0;  // for store num of vertex

    // there 4 type of render types indexed, triangle, fan, strip, for invalid case we make it defualt
    switch(type) {
        // if is indexed reander type
        case render_type::indexed:
        {
            // fill every vertex/indexed
            for (int i = 0; i < state.num_triangles; i++) {
                for (int j = 0; j < 3; j++) {
                    // fill the geo data for the indexed, increse by floats per vertex
                    input_geo[j].data = state.vertex_data + state.index_data[vertex_index] * state.floats_per_vertex;
                    vertex_index++; // next vertex
                    temporary_geo[j].data = input_geo[j].data; // store it to temp geo
                    state.vertex_shader(input_geo[j], temporary_geo[j], state.uniform_data); // then calculate the geo pos
                    tri_output[j] = temporary_geo[j];  // store our geo pos to our ouput
                }
                clip_triangle(state, tri_output[0], tri_output[1], tri_output[2], 0);
            }
            break;
        }
        // if is traingle reander type
        case render_type::triangle:
        {
            // calculate number of triangle, eac_edgeh traingle has 3 vertex
            num_triangles = state.num_vertices / 3;
            // fill eac_edgeh trangles
            for (int i = 0; i < num_triangles; i++) {
                for (int j = 0; j < 3; j++) { // fill eac_edgeh vertex data
                    // fill the geto datat of the traingle
                    input_geo[j].data = state.vertex_data + vertex_index; // triangle: increse each vertex
                    vertex_index += state.floats_per_vertex; // index of vetex increse
                    // calculate the geo data location
                    temporary_geo[j].data = input_geo[j].data;  // store it to temp geo
                    state.vertex_shader(input_geo[j], temporary_geo[j], state.uniform_data); // then calculate the geo pos
                    tri_output[j] = temporary_geo[j];  // store our data to ouput, for the 3 vector that send to clip traingle
                }
                // clip the traingle by call the function: internal state , 3 geo vector, and face start at 0
                clip_triangle(state, tri_output[0], tri_output[1], tri_output[2], 0);
            } 
            break;
        }
        // if is fan reander type
        case render_type::fan:
        {
            // number of tirangle in fan are vertices - 2
            num_triangles = state.num_vertices-2;
            for (int i = 0; i < num_triangles; i++) {
                for (int j = 0; j < 3; j++) {
                    // fill the geo data for the fan
                    if(j == 0) { // first vertex of the triangles
                        input_geo[j].data = state.vertex_data;  // we made it equal to first vertex, b/c every other vertex are connect to the first vertex
                    }
                    else { // if is 2, 3rd vertex of the traingle we do
                        input_geo[j].data = state.vertex_data + ( (state.floats_per_vertex) * (i + j) );
                    }
                    temporary_geo[j].data = input_geo[j].data; // store it to temp geo
                    state.vertex_shader(input_geo[j], temporary_geo[j], state.uniform_data); // then calculate the geo pos
                    tri_output[j] = temporary_geo[j];  // store our geo pos to our ouput
                }
                clip_triangle(state, tri_output[0], tri_output[1], tri_output[2], 0);
            }
            break;
        }
        // if is strip reander type
        case render_type::strip:
        {
            // for the strip the number of num_triangles is number of verices -2
            // for example: the strip that had 5 vertices, it has 3 num_triangles
            // hence, vertices = number of  num_triangles + 2 
            num_triangles = state.num_vertices-2;
            for (int i = 0; i < num_triangles; i++) {
                for (int j = 0; j < 3; j++) {
                    // fill the geo data for the strip
                    // depend on number of triangle and vertex of each traingle 
                    input_geo[j].data = state.vertex_data + ( (state.floats_per_vertex) * (i + j) ); // same as: &state.vertex_data[state.floats_per_vertex*(i+j)]; 
                    temporary_geo[j].data = input_geo[j].data; // store it to temp geo
                    state.vertex_shader(input_geo[j], temporary_geo[j], state.uniform_data); // then calculate the geo pos
                    tri_output[j] = temporary_geo[j];  // store our geo pos to our ouput

                }
                clip_triangle(state, tri_output[0], tri_output[1], tri_output[2], 0);
            }
            break;
        }
        // invaild type is our default case, since it shows invalid type is not the member of the render_type whe nwe run the program
        default:
        {
            std::cerr << "ERROR: Invalid render_type !!!" << std::endl;
            break;
        }
        
    }
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for eac_edgeh clipping face (face=0, 1, ..., 5) to
// clip against eac_edgeh of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    data_geometry temporary_geo[3] = { v0, v1, v2 }; // temporary store geo v0, v1, v2
    data_geometry a_vertex, b_vertex , c_vertex;  // for store the vertices
    data_geometry output_geo[3];  // for store our triangle to be clip

    // NDC: sign [1,-1]
    int sign = 0;  // 1 indicate inside, -1 indicate outside
    int num_vert = 0; // for store number of vertices that inside the image/plane
    // face clip 6 just rasterize the traingle
    bool is_inside[3] = {0}; // for store the bool number of the vertex a, b, c are inside(1) or not inside(0) 
    int gl_index = 0; // for store the index that makesure the gl position of x, y, z

    // if face is 6, we just rasterize the triangle
    if(face == 6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    //std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    // other face = 1,2,3,4,5, sign[1,-1], define there sign, purpose for finding there are outside or inside the image 
    // right, left plane against x
    // top, bottom plane against y
    // far, near plance against z
    // define x, y ,z(x = gl_Position[0],y = gl_Position[1], z =gl_Position[2])
    // face 0,1 for x, face 2, 3 for y, face 4,5 for z 
    // gl_Position[index], index = 0, 1, 2
    else if (face == 0) {  // right plane
        sign = 1;  // inside;
        gl_index = 0; // face 0 is for x, gl_position at index 0
    }
    else if (face == 1) {  // left plane
        sign = -1;  // ouside
        gl_index = 0; // face 0 is for -x, gl_position at index 0
    }
    else if (face == 2) {  // top plane
        sign = 1;  // insdie
        gl_index = 1; // face 2 is for y, gl_position at index 1
    }
    else if (face == 3) {  // bottom plane
        sign = -1;  // outside
        gl_index = 1; // face 3 is for -y, gl_position at index 1
    }
    else if (face == 4) { // far plane
        sign = 1;  // inside
        gl_index = 2; // face 4 is for z, gl_position at index 2
    }
    else if (face == 5) { // near plane
        sign = -1; // outside
        gl_index = 2; // face 5 is for -z, gl_position at index 2
    }

    // After defince each face's sign, we need to define is actually true inside and # of vertex
    // calculate the number of vertices are inside our image/plance
    for (int i = 0; i < 3; i++) {  // loop thur the 3 verteics 
        // deterimate # of vertices is inside the plane(note: sign=[1,-1])
        if (sign == 1) {  
            if (temporary_geo[i].gl_Position[gl_index] <= temporary_geo[i].gl_Position[3]) {  // if x(y,z) <= w then is inside, (x = gl_Position[0],y = gl_Position[1], z =gl_Position[2], w = gl_Position[3])
                is_inside[i] = 1;  // set our sign to 1 to indicate current vertex is inside
                num_vert++;  // number of vertex inside the image/plac_edgee add 1
            } 
        }
        else if (sign == -1) {  
            if (temporary_geo[i].gl_Position[gl_index] >= -temporary_geo[i].gl_Position[3]) { // if x(or y, or z) >= -w  then is inside
                is_inside[i] = 1;  // set our sign to 1 to indicate current vertex is inside
                num_vert++;  // number of vertex inside the image/plac_edgee add 1
            }
        }
    }

    // 4 cases: all vertices inside the image, all vertices outside the image, 1 vertex inside the image, 2 vertex inside the image
    // If all vertices are not in our image
    if (num_vert == 0) { // if number of vertex is zero, which means the triangle is completely outsdie of the image/plane
        return; //do nothing
    }
    // if all 3 vertices are inside
    else if (num_vert == 3) {
        clip_triangle(state, v0, v1, v2,face+1); // just clip the triangle, not extra jobs need to be done
        return;
    }
    // else if 1 or 2 vetrices insdie the plane, we need to do more
    // since we done for check for all vertices are inside and outside the image/plane
    // Below are for 2 more cases: only 1 vertex inside, and 2 vertices inside (we need to find the point along the edge of vertex that outside the image)
    for (int index = 0; index < 3; index++) {
        // (v0,v1,v2)
        // (v0,v0,v1) 
        // (v0,v1,v1)
        // if there are two 2 vertex and current index is outside the place(we need to find out the vertex that not inside the image)
        if ((num_vert == 2 && is_inside[index] == 0)) {
            // store the geo vector to vertex variable for calculate the point along the edge(ab,ac)
            a_vertex = temporary_geo[index];       // only store v0(v1,or v2) position, possible index (0,1,2)
            b_vertex = temporary_geo[(index+1)%3]; // only store v0 or v1 position, possible index(0,0,1)
            c_vertex = temporary_geo[(index+2)%3]; // only store v0 or v1 position, possible index(0,1,1)
        }
        // if only 1 vertex inside(there is one vertex and inside the plance)(we need to find out the vertex that insdie the image)
        else if ( (num_vert == 1 && is_inside[index] == 1)) {
             // store the get vertex to vertex variable for calculate the point along the edge(ab,ac) 
            a_vertex = temporary_geo[index];       // only store v0(v1,or v2) position, possible index(0,1,2)
            b_vertex = temporary_geo[(index+1)%3]; // only store v0 and v1 position, possible index(0,0,1)
            c_vertex = temporary_geo[(index+2)%3]; // only store v0 and v1 position , possible index(0,1,1)
        } 
    }

    // store the gl_Positon for calculate: alpha, beta, gamma 
    float a_position_W = a_vertex.gl_Position[3];
    float b_position_W = b_vertex.gl_Position[3];
    float c_position_W = c_vertex.gl_Position[3];
    // w depend on the sign to be negtive or positive, for e.g near sign = 1, w is pos, far sign = -1, w is neg, so on for the 4 more cases(left, right, top, and bottom plane)
    // find the β, γ, b/c we need to find the edge of ab(β), ac(γ)
    // e.g  β = (-w_b − x_b)/(x_a + w_a − w_b − x_b)    against left plane, sign = -1
    //      γ = (-w_c − x_c)/(x_a + w_a − w_c − x_c)    against bottom plane, sign = -1
    //      β = (w_b − x_b)/(x_a - w_a + w_b − x_b)     against right plane, sign = 1
    //      γ = (w_c − x_c)/(x_a - w_a + w_c − x_c)     against top plane, sign = 1
    float beta = ((sign * b_position_W) - b_vertex.gl_Position[gl_index]) / ((a_vertex.gl_Position[gl_index] - (sign * a_position_W)) + ((sign * b_position_W) - b_vertex.gl_Position[gl_index]));
    float gamma = ((sign * c_position_W) - c_vertex.gl_Position[gl_index]) / ((a_vertex.gl_Position[gl_index] - (sign * a_position_W)) + ((sign * c_position_W) - c_vertex.gl_Position[gl_index]));
    

    data_geometry ab_edge; // for points along the edge ab
    data_geometry ac_edge; // for points along the edge ac
    ab_edge.data = new float[state.floats_per_vertex];  // initlize the value
    ac_edge.data = new float[state.floats_per_vertex];  // initlize the value
    //  Implement color interpolation
    for(int i = 0; i < state.floats_per_vertex; i++)
    {
        // checking interp rules
        switch(state.interp_rules[i])
        {
            // if is invalid interp type display the error
            case interp_type::invalid:
            {
                std::cerr << "ERROR: Invalid interp_type!!!\n";
                break; 
            }
            // check interpolation type and interplotate the data depend on the type
            case interp_type::flat: // if is interpolation type of flat
            {
                // interpolate the data
                // store the all the data from the first vertex data
                ab_edge.data[i] = a_vertex.data[i];
                ac_edge.data[i] = a_vertex.data[i];
                break;
            }
            case interp_type::smooth: // if is interpolation type of smooth
            {
                // interpolate the data, correct the interpolation based on perspective-correct
                ab_edge.data[i]= beta * a_vertex.data[i] + ((1-beta)*b_vertex.data[i]);
                ac_edge.data[i]= gamma * a_vertex.data[i] + ((1-gamma)*c_vertex.data[i]);
                break;
            }
            case interp_type::noperspective: // if is interpolation type of noperspective
            {
                // interpolate the data, based on the image space
                float  k_valueue = (beta * a_position_W) + ((1-beta) * b_position_W);
                float alpha_temp = (beta * a_position_W)/k_valueue;
                ab_edge.data[i]= (alpha_temp * a_vertex.data[i]) + ((1-alpha_temp) * b_vertex.data[i]);
                float k_temp = (gamma * a_position_W) + ((1-gamma) * c_position_W);
                float gamma_temp = (gamma * a_position_W)/k_temp;
                ac_edge.data[i]= (gamma_temp * a_vertex.data[i]) + ((1-gamma_temp) * c_vertex.data[i]);
                break;
            }
            default: // invaild case is our default
            {
                std::cerr << "ERROR: Invalid interp_type!!!\n";
                break; 
            }   
        }
    }

    // store world bary coords for edge ab and ac
    for (int i = 0; i < 4; i++)
    {
        ab_edge.gl_Position[i]= (beta * a_vertex.gl_Position[i]) + ((1-beta)*b_vertex.gl_Position[i]);
        ac_edge.gl_Position[i]= (gamma * a_vertex.gl_Position[i]) + ((1-gamma)*c_vertex.gl_Position[i]);
    }


    // clip the triangles that only has one vertex insdie the image/plane
    if (num_vert == 1)
    {
        // one vertex and 2 vertics along with the edges(ab, ac)
        // traingle that have 3 coords, one completely inside, 2 coord along on the edge
        // traingle coords:(a coords, ac coords, ab coords)
        output_geo[0]=a_vertex;
        output_geo[1]=ab_edge;
        output_geo[2]=ac_edge;
        // just clip once at the time
        clip_triangle(state,output_geo[0],output_geo[1],output_geo[2],face+1);;
        return;
    }
    // clip the triangles that has 2 vertices insdie the image/plane
    else if (num_vert == 2)
    {
        // since two vertice inside, so we have 4 coords and we can divide to 2 traingles
        // so we need to clip 2 traingles when 2 vecrtices is inside
        // 2 vertices and on point along the edge
        // we have 2 case: 2 sub-traingles
        // sub-traingles that two vertices completely insdie and one vertex along on edge
        // traingles coords(ab coords, b coords, c coords)
        output_geo[0]=ab_edge;
        output_geo[1]=b_vertex;
        output_geo[2]=c_vertex;
        clip_triangle(state,output_geo[0],output_geo[1],output_geo[2],face+1);
        // sub-traingles that one vertex inside inside and two vertices along on edge
        // traingles second triangles coords(c coords, ac coords, ab coords)
        output_geo[0]=c_vertex;
        output_geo[1]=ac_edge;
        output_geo[2]=ab_edge;
        clip_triangle(state,output_geo[0],output_geo[1],output_geo[2],face+1);
        return;
    }
    /*else {
        //std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;

        // If all vertex in our image
        // which means the triangle is completely ousdie of the image/plane
        // x < w is outside
        if(a_vertex[2] < -a_vertex[3] && b_vertex[2] < -b_vertex[3] && c_vertex[2] < -c_vertex[3]) {
            return; // do nothing
        }
        // if there is one vertex is inside of the image/plane
        // a, b, or c is inside
        if ( a_vertex[2] >= -a_vertex[3] && b_vertex[2] < -b_vertex[3] && c_vertex[2] < -c_vertex[3] ) { // if only a is inside
            create_triangle2(state);
        }
        else if ( a_vertex[2] < -a_vertex[3] && b_vertex[2] >= -b_vertex[3] && c_vertex[2] < -c_vertex[3] ) { // if only b is inside

        }
        else if ( a_vertex[2] < -a_vertex[3] && b_vertex[2] < -b_vertex[3] && c_vertex[2] >= -c_vertex[3] ) { // if only c is inside

        }
        // if there are two vertices are inside of the image/pance
        // a and b; or a and c; or b and c;
        else if ( a_vertex[2] >= -a_vertex[3] && b_vertex[2] >= -b_vertex[3] && c_vertex[2] < -c_vertex[3] ) { // if a and b are inside

        }
        else if ( a_vertex[2] >= -a_vertex[3] && b_vertex[2] < -b_vertex[3] && c_vertex[2] >= -c_vertex[3] ) { // if a and c are inside

        }
        else if ( a_vertex[2] < -a_vertex[3] && b_vertex[2] >= -b_vertex[3] && c_vertex[2] >= -c_vertex[3] ) {  // if b and c is inside

        }
        // if all 3 vertices are inside
        else if(a_vertex[2] >= -a_vertex[3] && b_vertex[2] >= -b_vertex[3] && c_vertex[2] >= -c_vertex[3])
        {
            clip_triangle(state,v0,v1,v2,face+1);
        }
    }*/
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    //std::cout<<"TODO: implement rasterization"<<std::endl;
    // temporary store the 3 data_geomerty vector
    data_geometry temporary_geo_v[3] = { v0, v1, v2 };
    // constant max value for color
    const int MAX_COLOR_NUM = 255;
    // For fragment of the shading, allocate the array to hold the interpolation data
    float* frag_data = new float[MAX_FLOATS_PER_VERTEX];
    data_fragment fragment_data{frag_data}; 
    data_output frag_output;

    float x_coords[3];  // x is correspond to the x pixel coordinates for eac_edgeh vertex in the traingle
    float y_coords[3];  // y is correspond to the y pixel coordinates for eac_edgeh vertex in the traingle
    float z_coords[3];  // z hold the z-coordinates of the perspective transformation of eac_edgeh vertex in the triangle
    // for store the rangle of x and y pixel coordinates (range to be draw pixel with color)
    float x_min = 0;
    float x_max = 0;
    float y_min = 0;
    float y_max = 0;
    // hold the area of tri
    float area_ABC =0; 
    float area_PBC = 0;
    float area_APC = 0;
    float area_ABP = 0;
    // for hold the true alpha, gamma, beta
    float alpha = 0;
    float beta = 0;
    float gamma = 0; 
    // for store the value of half width and height of the image
    float width_half  = 0;
    float height_half = 0;
    // for store the coef to correct the interpolation
    float k_value = 0;

    // holde the W coord for our 3 geo vector
    float w_v0 = v0.gl_Position[3];
    float w_v1 = v1.gl_Position[3];
    float w_v2 = v2.gl_Position[3];

    // calculate the half width and half height of the image for calculate he x, y coords purpose
    width_half = state.image_width/2.0f;
    height_half = state.image_height/2.0f;

    // store our  x y z pos to our temp coordinate (NDC)
    // calcualte the pixel coordinates 
    for(int index = 0; index < 3; index++) 
    {
        x_coords[index] = (width_half  * (temporary_geo_v[index].gl_Position[0]/temporary_geo_v[index].gl_Position[3]) + (width_half  - 0.5f));
        y_coords[index] = (height_half * (temporary_geo_v[index].gl_Position[1]/temporary_geo_v[index].gl_Position[3]) + (height_half - 0.5f));
        z_coords[index] = temporary_geo_v[index].gl_Position[2]/temporary_geo_v[index].gl_Position[3];
    }

    // For draw the traingle, first we find the barcentric weights
    // below is calculate the total area of the triangle.
    //AREA(ABC) = 0.5((BxCy − CxBy) + (CxAy − AxCy) + (AxBy − BxAy))
    area_ABC = (0.5f * ((x_coords[1]*y_coords[2] - x_coords[2]*y_coords[1]) - (x_coords[0]*y_coords[2] - x_coords[2]*y_coords[0]) + (x_coords[0]*y_coords[1] - x_coords[1]*y_coords[0])));

    // calculate the min and max value our the x, y postion 
    // compare each coordinatie in x to find min x
    x_min = std::min(std::min(x_coords[0], x_coords[1]), x_coords[2]);
    // compare each coordinatie in x to find max x 
    x_max = std::max(std::max(x_coords[0], x_coords[1]), x_coords[2]);
    // compare each coordinatie in y to find min y
    y_min = std::min(std::min(y_coords[0], y_coords[1]), y_coords[2]);
    // compare each coordinatie in y to find max y
    y_max = std::max(std::max(y_coords[0], y_coords[1]), y_coords[2]);


    // for lecture ptt, there have algorithm for it
    // iterate through each pixel in range of x,y
    for(int j = y_min + 1; j < y_max; j++)
    {
        for(int i = x_min + 1; i < x_max; i++) 
        {
            // calcualte area of num_triangles 
            // AREA(PBC)= 0.5((BxCy − CxBy) + (By − Cy) ∗ i + (Cx − Bx) ∗ j)
            area_PBC = 0.5f * ((x_coords[1]*y_coords[2] - x_coords[2]*y_coords[1]) + (y_coords[1] - y_coords[2])*i + (x_coords[2] - x_coords[1])*j);
            // AREA(APC) = 0.5((CxAy − AxCy) + (Cy − Ay) ∗ i + (Ax − Cx) ∗ j)
            area_APC = 0.5f * ((x_coords[2]*y_coords[0] - x_coords[0]*y_coords[2]) + (y_coords[2] - y_coords[0])*i + (x_coords[0] - x_coords[2])*j);
            //AREA(ABP) = 0.5((AxBy − BxAy) + (Ay − By) ∗ i + (Bx − Ax) ∗ j)
            area_ABP = 0.5f * ((x_coords[0]*y_coords[1] - x_coords[1]*y_coords[0]) + (y_coords[0] - y_coords[1])*i + (x_coords[1] - x_coords[0])*j);
            // compute alpha, beta, gamma for (x,y)
            alpha = (area_PBC) / area_ABC;
            beta =  (area_APC) / area_ABC;
            gamma = (area_ABP) / area_ABC;

            // If the alpha , beta, gamma all in the triangle (0<=alpha_temp , beta, gamma<=1)
            if (alpha >= 0 && beta >= 0 && gamma >= 0) 
            { 
                // calculate the color, which calculate the z coordinate
                float temporary_geo_z_depth = alpha * z_coords[0] + beta * z_coords[1] + gamma * z_coords[2];
                // i + j + image_width is our pixel index
                unsigned pixel_index = i + j * state.image_width;  
                // Draw our pixel 
                // if the pixel is indide the traingle (pixel index = i + j * state.image_width)
                if(temporary_geo_z_depth < state.image_depth[pixel_index])
                {
                    //  for each float per vertex, we have to interpolate the data that 
                    //  depend on the interp_rule
                    for(int i = 0; i < state.floats_per_vertex; i++)
                    {
                        // interpolate the data that depend on the interpolate rule
                        switch(state.interp_rules[i])
                        {
                            // if is invalid interp type display the error
                            case interp_type::invalid:
                            {
                                std::cerr << "ERROR: Invalid interp_type!!!\n";
                                break; 
                            }
                            // for interpolation rule that is flat
                            case interp_type::flat:
                            {
                                // set all data to the first vertex data
                                fragment_data.data[i] = v0.data[i];
                                break;
                            }
                            // for the interpolation rule that is smooth
                            case interp_type::smooth:
                            {
                                // correct the interploation, based on the world space edge(perspective-correct)
                                // find the k value
                                // K = α/w_a+ β/w_b+ γ/w_c
                                k_value = (alpha / w_v0) + (beta  / w_v1) + (gamma / w_v2);
                                // correct the interpolation by correct the alpha, beta and gamma
                                float correct_alpha = alpha / (k_value * w_v0);  // alpha = (α/w_a)/k = α/(w_a*k)
                                float correct_beta  = beta  / (k_value * w_v1);  // beta = (β/w_b)/k = β/(w_a*k)
                                float correct_gamma = gamma / (k_value * w_v2);  // gamma = (γ/w_c)/k = β/(w_a*k)

                                // then interpolate the fragemt
                                fragment_data.data[i] = correct_alpha * v0.data[i] + correct_beta  * v1.data[i] + correct_gamma * v2.data[i];
                                break;
                            }
                            // for the interpolation that is noperspecitive
                            case interp_type::noperspective:
                            {
                                // interplation the fragment data depend on the (x,y)/image space
                                fragment_data.data[i] = alpha * v0.data[i] + beta  * v1.data[i] + gamma * v2.data[i];
                                break;
                            }
                            default: // invaild case is our default
                            {
                                std::cerr << "ERROR: Invalid interp_type!!!\n";
                                break; 
                            }        
                        }
                    }
                    // after interpolation data
                    // set the image depth to be z-buff coordinates depth
                    state.image_depth[pixel_index] = temporary_geo_z_depth;
                    // call the fragment shader with our interpolated data
                    state.fragment_shader(fragment_data, frag_output, state.uniform_data);
                    // then color it
                    state.image_color[pixel_index] = make_pixel((frag_output.output_color[0] * MAX_COLOR_NUM),(frag_output.output_color[1] * MAX_COLOR_NUM),(frag_output.output_color[2] * MAX_COLOR_NUM));
                }
            }
        }
    }
    // since we use allocate array, we need to delete the allocted memory and free the space.
    delete [] frag_data;
}
