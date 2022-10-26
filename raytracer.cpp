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
};

Vec3f operator*( const Vec3f& a,  const Vec3f& b)
{
    Vec3f result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

Vec3f operator*( const Vec3f& a, float b)
{
    Vec3f res;
    res.x = a.x * b;
    res.y = a.y * b;
    res.z = a.z * b;
    return res;
}

float dotProduct(const Vec3f &a, const Vec3f &b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}
Vec3f operator+( const Vec3f& a,  const Vec3f& b)
{
    Vec3f res;
    res.x = a.x + b.x;
    res.y = a.y + b.y;
    res.z = a.z + b.z;
    return res;
}

Vec3f operator-( const Vec3f& a,  const Vec3f& b)
{
    Vec3f res;
    res.x = a.x - b.x;
    res.y = a.y - b.y;
    res.z = a.z - b.z;
    return res;
}

Vec3f operator/( const Vec3f& a, float b)
{
    Vec3f res;
    res.x = a.x / b;
    res.y = a.y / b;
    res.z = a.z / b;
    return res;
}

Vec3f normalizeVector(const Vec3f& v)
{
    Vec3f res;
    float length = sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
    res.x = v.x / length;
    res.y = v.y / length;
    res.z = v.z / length;
    return res;
}

Vec3f hitPoint(const Ray& r, float t) //directly return res?
{

    
    return r.o + r.d * t;
}

Vec3f sphereNormal(const Vec3f& c, const Vec3f& p, const float& r)
{
    Vec3f res;
    res = normalizeVector( (p-c) / r );
    return res;
}

Vec3f triangleNormal(const Vec3f& v0, const Vec3f& v1, const Vec3f& v2)
{
    Vec3f res;
    res = (v1-v0)*(v2-v0);
    res = normalizeVector(res);
    return res;
}


Ray computeRay(const Camera& cam, int pixel_x, int pixel_y, int width, int height )
{

    Ray result;

    float l = cam.near_plane.x;
    float r = cam.near_plane.y;
    float b = cam.near_plane.z;
    float t = cam.near_plane.w;

    Vec3f v = cam.up;
    v = normalizeVector(v);
    Vec3f gaze = cam.gaze;
    gaze = normalizeVector(gaze);
    float dist = cam.near_distance;

    Vec3f u = gaze*v;
    u = normalizeVector(u);

    Vec3f m = cam.position + gaze*dist;

    Vec3f q = m + u*l + v*t;

    float su = (pixel_y+0.5)*(r-l)/width;
    float sv = (pixel_x+0.5)*(t-b)/height;

    Vec3f s = q+u*su-v*sv;

    result.o = cam.position;
    result.d = normalizeVector(s-cam.position);

    return result;

}

int discretizeColor(float color)
{
    int res;
    if(color > 255)
    {
        res = 255;
    }
    else
    {
        res = round(color) ;
    }
    return res;

}

float determinant(const Vec3f& v0, const Vec3f& v1, const Vec3f& v2 )
{
    return v0.x*(v1.y*v2.z - v1.z*v2.y) + v0.y*(v1.z*v2.x - v1.x*v2.z) + v0.z*(v1.x*v2.y - v1.y*v2.x);
}

float findDistance(const Vec3f& a, const Vec3f& b)
{
    return sqrtf(pow(a.x-b.x,2) + pow(a.y-b.y,2) + pow(a.z-b.z,2));
}

float findLength(const Vec3f& a)
{
    return sqrtf(a.x*a.x + a.y*a.y + a.z*a.z);
}

void checkSphereIntersection(const Ray& ray,  const Scene& scene, const int& i, float& t_min, int& type_of_closest_object, int& index_of_closest_object )
{
    Vec3f d = ray.d;
    Vec3f o = ray.o;

    Sphere sphere = scene.spheres[i];

    Vec3f c = scene.vertex_data[sphere.center_vertex_id-1];

    float r = sphere.radius;

    float A = dotProduct(ray.d, ray.d);
	Vec3f tmp = ray.o - c;
	float B = 2*dotProduct(ray.d, tmp);
	float C = dotProduct(tmp, tmp)-r*r;

	float disc = B*B - 4*A*C;

    if(disc<0)
    {
        return;
    }
    else // disc = 1 & 2. use epsilon if case 1 to be implemented individually.
    {
        float t1 = (-1*B+sqrt(disc))/(2*A);
        float t2 = (-1*B-sqrt(disc))/(2*A);
        float tmin = fmin(t1,t2);


        if ( tmin < t_min && tmin > 0  )
        {

            t_min = tmin;
            type_of_closest_object = 0;
            index_of_closest_object = i;
        }

        return;
    }

}

void checkTriangleIntersection(const Ray& ray,  const Scene& scene, const int& i, float& t_min, int& type_of_closest_object, int& index_of_closest_object )
{
    Vec3f d = ray.d;
    Vec3f o = ray.o;

    Triangle triangle = scene.triangles[i];
    Vec3f v0 = scene.vertex_data[triangle.indices.v0_id -1 ];
    Vec3f v1 = scene.vertex_data[triangle.indices.v1_id -1 ];
    Vec3f v2 = scene.vertex_data[triangle.indices.v2_id -1 ];

    float detA,beta,gamma,t;

    detA=determinant(v0-v1, v0-v2, d);
    beta = determinant(v0-o,v0-v2,d)/detA;
    gamma = determinant(v0-v1,v0-o,d)/detA;
    t = determinant(v0-v1,v0-v2,v0-o)/detA;

    if(beta+gamma > 1)
    {
        return;
    }
    else if(beta<0)
    {
        return;
    }
    else if(gamma<0)
    {
        return;
    }
    else //Ray hit with some triangle. Check if it is the closest.
    {
        if(t<t_min && t > 0)
        {
            t_min=t;
            type_of_closest_object=1;
            index_of_closest_object=i;
            return;
        }
        else
        {
            return;
        }
    }
}

bool hitBoundingBox(const Ray& ray, const Vec3f& c, float& r)
{

    Vec3f d = ray.d;
    Vec3f o = ray.o;

    float A = dotProduct(ray.d, ray.d);
	Vec3f tmp = ray.o - c;
	float B = 2*dotProduct(ray.d, tmp);
	float C = dotProduct(tmp, tmp)-r*r;

	float disc = B*B - 4*A*C;

    if(disc<0)
    {
        return false;
    }
    return true;
}


void checkMeshIntersection(const Ray& ray,  const Scene& scene, const int& i, int& face_number, float& t_min, int& type_of_closest_object, int& index_of_closest_object )
{
    Vec3f d = ray.d;
    Vec3f o = ray.o;
    
    Mesh mesh = scene.meshes[i];
    int face_amount = mesh.faces.size();
    /*
    float x_min, x_max, y_min, y_max, z_min, z_max;
    x_min = y_min = z_min = INFTY;
    x_max = y_max = z_max = 0;
    
    for(int j=0; j<face_amount;j++)
    {
        
        Vec3f v0 = scene.vertex_data[ mesh.faces[j].v0_id -1 ];
        Vec3f v1 = scene.vertex_data[ mesh.faces[j].v1_id -1 ];
        Vec3f v2 = scene.vertex_data[ mesh.faces[j].v2_id -1 ];
        float xmin, xmax, ymin, ymax, zmin, zmax;
        xmin = min(v0.x,min(v1.x,v2.x));
        if(xmin < x_min) x_min=xmin;

        xmax = max(v0.x,max(v1.x,v2.x));
        if(xmax > x_max) x_max = xmax;

        ymin = min(v0.y,min(v1.y,v2.y));
        if(ymin < y_min) y_min=ymin;

        ymax = max(v0.y,max(v1.y,v2.y));
        if(ymax > y_max) y_max = ymax;

        zmin = min(v0.z,min(v1.z,v2.z));
        if(zmin < z_min) z_min=zmin;

        zmax = max(v0.z,max(v1.z,v2.z));
        if(zmax > z_max) z_max = zmax;
        
    }
    
    Vec3f c;
    float r;

    c.x = (x_min + x_max)/2.0;
    c.y = (y_min + y_max)/2.0;
    c.z = (z_min + z_max)/2.0;
    
    r = max( (abs(x_min - x_max))/2.0, max( (abs(y_min - y_max))/2.0,  (abs(z_min - z_max))/2.0) );

    
    if (!hitBoundingBox(ray, c, r))
    {
        return;
    }
    */
    for(int j=0; j<face_amount; j++)
    {

        Vec3f v0 = scene.vertex_data[ mesh.faces[j].v0_id -1 ];
        Vec3f v1 = scene.vertex_data[ mesh.faces[j].v1_id -1 ];
        Vec3f v2 = scene.vertex_data[ mesh.faces[j].v2_id -1 ];
        float detA,beta,gamma,t;

        detA=determinant(v0-v1, v0-v2, d);
        beta = determinant(v0-o,v0-v2,d)/detA;
        gamma = determinant(v0-v1,v0-o,d)/detA;
        t = determinant(v0-v1,v0-v2,v0-o)/detA;

        if(beta+gamma > 1)
        {
            continue;
        }
        else if(beta<0)
        {
            continue;
        }
        else if(gamma<0)
        {
            continue;
        }
        else //Ray hit with some triangle. Check if it is the closest.
        {
            if(t<t_min && t > 0)
            {
                t_min=t;
                type_of_closest_object=2;
                index_of_closest_object=i;
                face_number = j;
                continue;
            }
            else
            {
                continue;
            }
        }
    }

}


void getEffects( const Scene& scene, Vec3f& color, const Vec3f& hit_location, const Vec3f& camera_location, const Vec3f& normal, const Vec3f& light_location, const Vec3f& k_d, const Vec3f& k_s, const Vec3f& intensity, 
                float& phong_exponent, float& shadow_epsilon, int& sphere_count, int& triangle_count, int& mesh_count )
{
    
    Vec3f w_i, w_o, h, irradiance;
    float cos_theta, cos_alpha, squared_distance_to_light;

    w_i = normalizeVector(light_location - hit_location);

    /* Shadow Check */
    
    Ray shadowRay;
    shadowRay.o = normal*shadow_epsilon + hit_location;
    shadowRay.d = w_i;
    int x = 0;
    float t = INFTY;
    int type_of_closest_object = -1;
    int index_of_closest_object = -1;
    int face_num = -1;

    for( int a = 0; a < sphere_count ; a++) //for each sphere get t_min
    {
        checkSphereIntersection(shadowRay, scene, a , t,  type_of_closest_object, index_of_closest_object);
        if(t<INFTY)
        {
            return;
        }
    }

    for(int a = 0; a < triangle_count ; a++) //for each triangle get t_min
    {
        checkTriangleIntersection(shadowRay, scene, a, t, type_of_closest_object, index_of_closest_object);
        if(t<INFTY)
        {
            return;
        }
    }

    for(int a = 0; a< mesh_count ; a++) //for each mesh get t_min
    {
        checkMeshIntersection(shadowRay, scene, a, face_num, t, type_of_closest_object, index_of_closest_object);
        if(t<INFTY)
        {
            return;
        }
    }

    /* Shadow Check */


    cos_theta =  max(0.0f,dotProduct(w_i,normal));

    squared_distance_to_light = pow(findDistance(hit_location, light_location),2);

    w_o = normalizeVector(camera_location - hit_location);
    h = normalizeVector(w_i+w_o);
    cos_alpha = max(0.0f, dotProduct(normal,h));


    irradiance.x = intensity.x/squared_distance_to_light;
    irradiance.y = intensity.y/squared_distance_to_light;
    irradiance.z = intensity.z/squared_distance_to_light;

    color.x += k_d.x * cos_theta * irradiance.x;
    color.y += k_d.y * cos_theta * irradiance.y;
    color.z += k_d.z * cos_theta * irradiance.z;

    color.x += k_s.x * pow(cos_alpha,phong_exponent) * irradiance.x;
    color.y += k_s.y * pow(cos_alpha,phong_exponent) * irradiance.y;
    color.z += k_s.z * pow(cos_alpha,phong_exponent) * irradiance.z;
    
    
}

Vec3f getColor(const Scene& scene, const Vec3f& hit_location, const Vec3f& camera_location, const Vec3f& normal, int& index, int& type, int& sphere_count, int& triangle_count, int& mesh_count )
{
    
    Vec3f  light_location;
    int amount_of_point_lights = scene.point_lights.size();
    float shadow_ray_epsilon = scene.shadow_ray_epsilon;

    Vec3f color, w_i, w_o, h, k_d, k_s, irradiance, intensity, ambient;
    float cos_theta, cos_alpha, squared_distance_to_light, phong_exponent;

    switch (type)
    {
    case 0:/* sphere  */

        ambient = scene.materials[scene.spheres[index].material_id -1 ].ambient;

        color.x = scene.ambient_light.x  * ambient.x;
        color.y = scene.ambient_light.y  * ambient.y;
        color.z = scene.ambient_light.z  * ambient.z;

        for(int i=0; i < amount_of_point_lights; i++)
        {
            

            light_location = scene.point_lights[i].position;
            k_d = scene.materials[scene.spheres[index].material_id -1 ].diffuse;
            phong_exponent = scene.materials[scene.spheres[index].material_id -1 ].phong_exponent;
            k_s = scene.materials[scene.spheres[index].material_id -1 ].specular;
            intensity = scene.point_lights[i].intensity;

            getEffects(scene, color, hit_location, camera_location, normal, light_location, k_d, k_s, intensity, phong_exponent, shadow_ray_epsilon, sphere_count, triangle_count, mesh_count  );
       
        }
        
        break;

    case 1:/* triangle  */
        ambient = scene.materials[scene.triangles[index].material_id -1 ].ambient;

        color.x = scene.ambient_light.x  * ambient.x;
        color.y = scene.ambient_light.y  * ambient.y;
        color.z = scene.ambient_light.z  * ambient.z;

        for(int i=0; i < amount_of_point_lights; i++)
        {
            

            light_location = scene.point_lights[i].position;
            k_d = scene.materials[scene.triangles[index].material_id -1 ].diffuse;
            phong_exponent = scene.materials[scene.triangles[index].material_id -1 ].phong_exponent;
            k_s = scene.materials[scene.triangles[index].material_id -1 ].specular;
            intensity = scene.point_lights[i].intensity;

            getEffects(scene, color, hit_location, camera_location, normal, light_location, k_d, k_s, intensity, phong_exponent, shadow_ray_epsilon, sphere_count, triangle_count, mesh_count  );       
        

        }
        break;

    case 2:/* mesh  */
        ambient = scene.materials[scene.meshes[index].material_id -1 ].ambient;

        color.x = scene.ambient_light.x  * ambient.x;
        color.y = scene.ambient_light.y  * ambient.y;
        color.z = scene.ambient_light.z  * ambient.z;

        for(int i=0; i < amount_of_point_lights; i++)
        {
            

            light_location = scene.point_lights[i].position;
            k_d = scene.materials[scene.meshes[index].material_id -1 ].diffuse;
            phong_exponent = scene.materials[scene.meshes[index].material_id -1 ].phong_exponent;
            k_s = scene.materials[scene.meshes[index].material_id -1 ].specular;
            intensity = scene.point_lights[i].intensity;

            getEffects(scene, color, hit_location, camera_location, normal, light_location, k_d, k_s, intensity, phong_exponent, shadow_ray_epsilon, sphere_count, triangle_count, mesh_count  );
       
        
        }
        break;

    case -1:
        color.x = scene.background_color.x ;
        color.y = scene.background_color.y ;
        color.z = scene.background_color.z ;
        break;
    }
    return color;
}






int main(int argc, char* argv[])
{

    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    int camera_amount = scene.cameras.size();

    for(int i=0 ; i < camera_amount ; i++)
    {

        int width = scene.cameras[i].image_width;
        int height = scene.cameras[i].image_height;
        unsigned char* image = new unsigned char [width * height * 3];

        Camera c = scene.cameras[i];
        Vec3f camera_location = c.position;
        int sphere_count = scene.spheres.size();
        int triangle_count = scene.triangles.size();
        int mesh_count = scene.meshes.size();
        int ppm_pixel=0;
        for(int j=0 ; j < height ; j++)
        {
            for(int k=0 ; k < width ; k++)
            {
                int type_of_closest_object = -1;
                int index_of_closest_object = -1;
                int face_num = -1;
                Ray ray = computeRay(c, j, k, width, height);
                Vec3f final_color;
                float t_min = INFTY;
                for( int a = 0; a < sphere_count ; a++) //for each sphere get t_min
                {
                    checkSphereIntersection(ray, scene, a , t_min,  type_of_closest_object, index_of_closest_object);
                }

                for(int a = 0; a < triangle_count ; a++) //for each triangle get t_min
                {
                    checkTriangleIntersection(ray, scene, a, t_min, type_of_closest_object, index_of_closest_object);
                }

                for(int a = 0; a< mesh_count ; a++) //for each mesh get t_min
                {
                    checkMeshIntersection(ray, scene, a, face_num, t_min, type_of_closest_object, index_of_closest_object);
                }

                Vec3f hit_point = hitPoint(ray,t_min);
                Vec3f normal;
                if(type_of_closest_object == 0)
                {
                    Sphere s = scene.spheres[index_of_closest_object];
                    normal = sphereNormal(scene.vertex_data[s.center_vertex_id -1], hit_point, s.radius);
                    
                }
                else if(type_of_closest_object == 1)
                {
                    Triangle t = scene.triangles[index_of_closest_object];
                    normal = triangleNormal(scene.vertex_data[t.indices.v0_id-1], scene.vertex_data[t.indices.v1_id-1], scene.vertex_data[t.indices.v2_id-1]);
                   
                }
                else if(type_of_closest_object == 2)
                {
                    Face f = scene.meshes[index_of_closest_object].faces[face_num];
                    normal = triangleNormal(scene.vertex_data[f.v0_id-1], scene.vertex_data[f.v1_id-1], scene.vertex_data[f.v2_id-1]);
                    
                }
                
                final_color = getColor(scene, hit_point, camera_location, normal, index_of_closest_object,type_of_closest_object, sphere_count, triangle_count, mesh_count);
        
                
                image[ppm_pixel++] =  discretizeColor(final_color.x);
                image[ppm_pixel++] =  discretizeColor(final_color.y);
                image[ppm_pixel++] =  discretizeColor(final_color.z);

            }

        }
        write_ppm(c.image_name.c_str(), image, width, height);
    }

}
