#include <iostream>
#include "parser.h"
#include "ppm.h"

#include "math.h"
#define INFTY 99999

using namespace std;
using namespace parser;



typedef unsigned char RGB[3];

struct Ray
{
    Vec3f o;
    Vec3f d;
    int t;
};



Vec3f operator*(const Vec3f& a, const Vec3f& b)
{
    Vec3f result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

Vec3f operator*(const Vec3f& a, float b)
{
    Vec3f res;
    res.x = a.x * b;
    res.y = a.y * b;
    res.z = a.z * b;
    return res;
}

float dotProduct(Vec3f& a, Vec3f& b)
{
 
    float result = 0;
 
    result += (a.x * b.x) + (a.y * b.y) + (a.z + b.z); 

    return result;
}

Vec3f operator+(const Vec3f& a, const Vec3f& b)
{
    Vec3f res;
    res.x = a.x + b.x;
    res.y = a.y + b.y;
    res.z = a.z + b.z;
    return res;
}

Vec3f operator-(const Vec3f& a, const Vec3f& b)
{
    Vec3f res;
    res.x = a.x - b.x;
    res.y = a.y - b.y;
    res.z = a.z - b.z;
    return res;
}

Vec3f normalizeVector(Vec3f& v)
{
    Vec3f res;
    float length = sqrt( pow(v.x,2) + pow(v.y,2) + pow(v.z,2) );
    res.x = v.x / length;
    res.y = v.y / length;
    res.z = v.z / length;
    return res;
}


Ray computeRay(Camera& cam, int& pixel_x, int& pixel_y, int& width, int& height )
{
    Ray result;  
    
    float l = cam.near_plane.x;
    float r = cam.near_plane.y;
    float b = cam.near_plane.z;
    float t = cam.near_plane.w;
    float dist = cam.near_distance;


    Vec3f gaze = normalizeVector( cam.gaze );
    Vec3f v = cam.up;
    Vec3f w = gaze*(-1);
    Vec3f u = v*w;

    Vec3f manhattan_distance = gaze*dist;
    Vec3f m = cam.position + manhattan_distance;

    Vec3f q = m + u*l + v*t;
    
    float su = (pixel_x+0.5)*(r-l)/width;
    float sv = (pixel_y+0.5)*(t-b)/height;

    Vec3f s = q + u*su - v*sv ;

    Vec3f temp = s - cam.position;
    result.o = cam.position;
    result.d = normalizeVector(temp);
    t = INFTY;

    return result;
}

void checkIntersection(Ray& ray,  Scene& scene, int& type_of_object)
{
    
}


int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.
    //
    // Normally, you would be running your ray tracing
    // code here to produce the desired image.
    int camera_amount = scene.cameras.size();
    Vec3i background_color = scene.background_color;
    for(int i=0 ; i < camera_amount ; i++)
    {

        string output_file_name = scene.cameras[i].image_name;
        int width = scene.cameras[i].image_width;
        int height = scene.cameras[i].image_height;
        Camera c = scene.cameras[i];

        for(int j=0 ; j < width ; j++)
        {
            for(int k=0 ; k < height ; k++)
            {
                int t_min = INFTY;
                int mesh_count = scene.meshes.size();
                int sphere_count = scene.spheres.size();
                int triangle_count = scene.triangles.size();

                Ray ray = computeRay(c, j, k, width, height);
                
                for( int a = 0; a < sphere_count ; a++) //for each sphere get t_min
                {
                    checkIntersection(ray, scene, 0);
                }

                for(int a = 0; a < triangle_count ; a++) //for each triangle get t_min
                {

                }
                for(int a = 0; a< mesh_count ; a++) //for each mesh get t_min
                {

                }
                

            }

        }
    }



    const RGB BAR_COLOR[8] =
    {
        { 255, 255, 255 },  // 100% White
        { 255, 255,   0 },  // Yellow
        {   0, 255, 255 },  // Cyan
        {   0, 255,   0 },  // Green
        { 255,   0, 255 },  // Magenta
        { 255,   0,   0 },  // Red
        {   0,   0, 255 },  // Blue
        {   0,   0,   0 },  // Black
    };

    int width = 640, height = 480;
    int columnWidth = width / 8;

    unsigned char* image = new unsigned char [width * height * 3];

    int i = 0;
    for (int y = 0; y < height; ++y)
    {
        for (int x = 0; x < width; ++x)
        {
            int colIdx = x / columnWidth;
            image[i++] = BAR_COLOR[colIdx][0];
            image[i++] = BAR_COLOR[colIdx][1];
            image[i++] = BAR_COLOR[colIdx][2];
        }
    }

    write_ppm("test.ppm", image, width, height);

}
