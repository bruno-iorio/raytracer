#define _CRT_SECURE_NO_WARNINGS 1
#define RAND_MAX 1e20
#define PI  3.14159265358979
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION

#include "stb_image_write.h"
#include "stb_image.h"
#include <iostream>
#include <vector>
#include <ctime>
#include <omp.h>


double sqr(double x){return x*x;}
double mod(double x){return (x > 0) ? x : -x;}

class Vector {
public:
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    double norm2() const {
        return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
    }
    double norm() const {
        return sqrt(norm2());
    }
    void normalize() {
        double n = norm();
        data[0] /= n;
        data[1] /= n;
        data[2] /= n;
    }
    double operator[](int i) const { return data[i]; };
    double& operator[](int i) { return data[i]; };
    double data[3];
};
 
Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const double b) {
    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

class LightSource{
    public:
        LightSource(const Vector& pos, double I)
        : pos(pos), I(I) {}

        LightSource()
        : pos(Vector(0,0,0)), I(1e8) {}

        Vector pos;
        double I;
};

class Ray {
public:
    Ray(const Vector& O ,const Vector& u) 
    : O(O),u(u) {}
    Ray()
    : O(Vector(0,0,0)),u(Vector(0,0,0)){}

    Vector O;
    Vector u;
};
 


class Sphere {
public:
    Sphere(const Vector& C, double R,Vector albedo,bool mirror,bool transparent, double n)
    : C(C), R(R), albedo(albedo), mirror(mirror), transparent(transparent), n(n) {}

    Sphere()
    : C(Vector(0,0,0)), R(0), albedo(Vector(0,0,0)), mirror(false), transparent(false), n(1) {}

    bool intersect(const Ray& r, Vector& pos, bool in){
        double delta = sqr(dot(r.u,r.O - C)) - dot(r.O - C, r.O - C) + R*R;
        if (delta < 0) return false;
        double x = dot(r.u,C-r.O);
        double t1 = x - sqrt(delta);
        double t2 = x + sqrt(delta);
        if (t2 < 0) return false;
        if (!in) {pos = r.O + t1*r.u;}
        else {pos = r.O + t2*r.u;}
        return true;
    }

    Vector C;
    Vector albedo;
    double R, n;
    bool mirror, transparent;

};
 

class Scene{
    public: 
        Scene(const std::vector<Sphere>& objects, const LightSource& light_source,double n)
        : objects(objects), light_source(light_source), n(n){}

        void add_object(Sphere& obj){objects.push_back(obj);}

        bool detect_shadow(Vector& P1,Vector& N,Sphere& s1, Vector& P2,bool &in){
            Vector L = light_source.pos;
            Vector N2(0,0,0);
            Ray r(P1 + 1e-4*N,(L-P1)/(L-P1).norm());

            double d,c1,c2,c3;

            Sphere s(Vector(0,0,0),0,Vector(0,0,0),false,false,1);
            if (this->intersect(r,s,P2,in)){
                if(std::pow((P2 - P1).norm(),2) <= std::pow((L-P2).norm(),2)){
                    return true;
                }
                return false;
            }
            return false;
        }

    bool get_color(Ray r,int raydepth ,double& color1, double& color2, double& color3,bool &in){
            if(raydepth == 0){
                return false;
            }
            
            Vector P,N;
            Sphere s;
            
            bool res = this->intersect(r,s,P,in);
            double I = light_source.I;

            Vector L = light_source.pos;
            Vector ro = s.albedo;
            Vector shadow_inter;
            N = P - s.C;
            N.normalize();

            if(s.mirror){ //if mirror just computes reflection
                Ray r_ref(P + 1e-4* N, r.u - 2*dot(r.u,N) * N);
                res = this->get_color(r_ref,raydepth-1,color1,color2,color3,in);
                return res;
            }
            else if(s.transparent){ //annoyting as hell!
                                    // dot(N,r.u) < 0 := (N and r.u have different orientation) == r.u is entering the ball
              Vector w_i = r.u;
              w_i.normalize();
              // fresnel reflection things:
              double n1 = (dot(N, r.u)) < 0 ? this->n : s.n;
              double n2 = (dot(N, r.u)) < 0 ? s.n : this->n;
              double k_0 = sqr(n1-n2)/sqr(n1+n2);
              double R = k_0 + (1 - k_0)*(1-mod(dot(N,w_i)));
              double T = 1 - R;

              bool in1 = (dot(N,r.u) < 0) ? true : false;
              bool in2 = (dot(N,r.u) < 0) ? false : true;

              Vector N1(0,0,0), N2(0,0,0);
              Vector P1 = in1 ? P - 1e-4 * N : P + 1e-4 * N;
              Vector P2 = in2 ? P - 1e-4 * N : P + 1e-4 * N;

              if (1 - sqr(n1/n2)*(1 - sqr(dot(w_i,N))) < 0){ // total reflection occurs
                Ray r_reflec(P2, w_i - 2*dot(w_i,N) * N);
                return this->get_color(r_reflec, raydepth-1,color1,color2,color3,in2);
              }
              else{
                Vector w_t = (n1/n2) * (w_i - dot(w_i,N)*N);
                Vector w_n = N * (sqrt(1 - sqr(n1/n2)*(1 - sqr(dot(w_i,N)))));
                Vector w = in1 ? w_t - w_n : w_t + w_n;
                Ray r_refrac(P1, w);
                Ray r_reflec(P2 ,w_i - 2*dot(w_i,N)*N);

                if(((double) rand()) / ((double) RAND_MAX) <= T){ // refracted ray :
                  res &= this->get_color(r_refrac,raydepth-1,color1,color2,color3,in1);
                }
                else{ //reflected ray :
                  res &= this->get_color(r_reflec,raydepth-1,color1,color2,color3,in2);
                }
                return res;
              }
            
            }

            else{
              color1 = I/(4 * PI * dot(L-P,L-P)) * ro[0]/PI * dot(N,(L - P)/sqrt(dot(L - P,L - P)));
              color2 = I/(4 * PI * dot(L-P,L-P)) * ro[1]/PI * dot(N,(L - P)/sqrt(dot(L - P,L - P)));
              color3 = I/(4 * PI * dot(L-P,L-P)) * ro[2]/PI * dot(N,(L - P)/sqrt(dot(L - P,L - P)));
            }
            
            if(dot(N,(L-P)/sqrt(dot(L-P,L-P))) < 0 || this->detect_shadow(P,N,s,shadow_inter,in)){
                color1 = 0;
                color2 = 0;
                color3 = 0;
                res = false;
            }
        return res;
    }

    bool intersect(Ray r, Sphere& s, Vector& P, bool &in){
        double distance = 1e10;
        bool found = false;
        Vector P1, origin(0,0,0);
        std::vector<Sphere>::iterator id = objects.begin(); 

        for(std::vector<Sphere>::iterator i = objects.begin(); i != objects.end(); i++){
            Vector intersection(0,0,0);
            double d;
            if(i->intersect(r,intersection,in)){ //checks intersection for each sphere
                double d = (intersection - origin).norm(); // computes the distance
                if(distance){
                    if (d < distance){ // if this distance is smaller than the global min distance, than it becomes the smalles distance and we change the object id
                        distance = d;
                        id = i;
                        P1 = intersection;
                        found = true;
                    }
                }
            }
        }
        if (found){
            s = *id;
            P = P1;
        }
        return found;
    }

    double n; 
    std::vector<Sphere> objects;
    LightSource light_source;

};


int main() {
    srand(time(0));
    double gamma=2.2;
    int W = 512;
    int H = 512;
    int rays_pixel = 100;
    Vector Camera_origin(0,0,55);
    Vector origin(0,0,0);

    LightSource light_source(Vector(-10,20,40), 1.4*1e10);

    double fov = PI* 60.0/180.0;
    double d = -W/(2*tan(fov/2));

    // Creating objects
    Sphere main_sphere1(Vector(0,0,0),10,Vector(0.8,0.8,0.8),false,true,1.5);
    Sphere main_sphere2(Vector(21,0,0),10,Vector(0.8,0.8,0.8),true,false,1.33);
    Sphere main_sphere3(Vector(-21,0,0),10,Vector(0.8,0.8,0.8),false,false,1.5);
    Sphere main_sphere4(Vector(-21,0,0),0.7,Vector(0.8,0.8,0.8),false,false,1);
    Sphere wall_left(Vector(-1000,0,0),940,Vector(0.9,0.2,0.9),false,false,1);
    Sphere wall_front(Vector(0,0,1000),940,Vector(0.9,0.4,0.3),false,false,1);
    Sphere wall_back(Vector(0,0,-1000),940,Vector(0.4,0.8,0.7),false,false,1);
    Sphere floor(Vector(0,-1000,0),990,Vector(0.3,0.4,0.7),false,false,1);
    Sphere ceiling(Vector(0,1000,0),940,Vector(0.2,0.5,0.9),false,false,1);
    Sphere wall_right(Vector(1000,0,0),940,Vector(0.6,0.5,0.1),false,false,1);

    std::vector<Sphere> objects = {main_sphere1,main_sphere2,main_sphere3,wall_back,wall_left,wall_front,wall_right,ceiling,floor};

    std::vector<unsigned char> image(W * H * 3, 0);

    Scene scene(objects,light_source,1);
    omp_set_num_threads(rays_pixel);

    std::vector<unsigned char> image1[rays_pixel];
    for (int i = 0; i < rays_pixel; i++){
      image1[i] = std::vector<unsigned char>(W * H * 3, 0);
    }

    #pragma omp parallel
    {
    int thread_id = omp_get_thread_num();
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector dir_vec(j-W/2 + 0.5, H/2 - i + 0.5, d);
            dir_vec.normalize();
            Ray dir_ray(Camera_origin,dir_vec);
            double color1=0,color2=0,color3=0;
            bool cam_env = false;
            if (scene.get_color(dir_ray,100,color1,color2,color3,cam_env)){
              image1[thread_id][(i * W + j) * 3 + 0] = (unsigned char) std::min(255.0, std::pow(color1,1/gamma));
              image1[thread_id][(i * W + j) * 3 + 1] = (unsigned char) std::min(255.0, std::pow(color2,1/gamma));
              image1[thread_id][(i * W + j) * 3 + 2] = (unsigned char) std::min(255.0, std::pow(color3,1/gamma));
          }
        }
      }
    }

    for (int i = 0; i < H; i++){
      for (int j = 0; j < W; j++){
        int k = 0;
        double color1=0,color2=0,color3=0;
        for (int id = 0 ; id < rays_pixel; id++){
             color1 += image1[id][(i * W + j) * 3 + 0];
             color2 += image1[id][(i * W + j) * 3 + 1];
             color3 += image1[id][(i * W + j) * 3 + 2];
        }
        image[(i * W + j) * 3 + 0] = color1/(double)rays_pixel;
        image[(i * W + j) * 3 + 1] = color2/(double)rays_pixel;
        image[(i * W + j) * 3 + 2] = color3/(double)rays_pixel;
        }
      }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
    return 0;
}

