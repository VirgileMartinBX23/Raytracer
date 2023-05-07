#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <typeinfo> 

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "Vector.h"
#include "Ray.h"
#include "Geometry.h"
#include "Scene.h"

#define PI 3.141592653589793238

int main() {
    int W = 256;
    int H = 256;

    Vector Q(0,0,55); //camera center
    double alpha = 60. * PI / 180.; //fov of 60° in rad 
    double epsilon = 0.001; //for the numerical precision
    int depth = 5; //ray depth
    int repetition = 32; //paths per pixel

    Scene scene;
    scene.I = 2E10; 
    scene.S = Vector(-10, 20, 40);

    //walls 
    Sphere wall_l(Vector(-1000,0,0), 940, Vector(0.5, 0.8, 0.1), false, false, 1.);
    Sphere wall_r(Vector(1000,0,0), 940, Vector(0.9, 0.2, 0.3), false, false, 1.);
    Sphere wall_t(Vector(0,1000,0), 940, Vector(0.3, 0.5, 0.3), false, false, 1.);
    Sphere wall_b(Vector(0,-1000,0), 990, Vector(0.6, 0.5, 0.7), false, false, 1.);
    Sphere wall_f(Vector(0,0, -1000), 940, Vector(0.1, 0.6, 0.7), false, false, 1.);
    Sphere wall_be(Vector(0,0, 1000), 940, Vector(1., 0.75, 0.8), false, false, 1.);
    scene.add_sphere(&wall_l);
    scene.add_sphere(&wall_r);
    scene.add_sphere(&wall_t);
    scene.add_sphere(&wall_b);
    scene.add_sphere(&wall_f);
    scene.add_sphere(&wall_be);

    //objects 
    //Sphere ball(Vector(0,0,0), 10, Vector(1., 1., 1.), false, true, 1.4);

    TriangleMesh Cat(Vector(1.,1.,1.), false, false, 1.);
    Cat.readOBJ("cat.obj");
    Cat.translated_and_scaled(Vector(0,-10,0), 0.6);
    Cat.compute_the_bounding_box();

    //scene.add_sphere(&ball);
    scene.add_mesh(&Cat);
 
    std::vector<unsigned char> image(W * H * 3, 0);
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {
            Vector color;
            double n1 = 1. ;//refraction index of air 
            //compute the color of the pixel
            for (int k = 0; k < repetition; k++){
                //compute the ray direction from the camera fro each pixel
                Vector ray_u(j - W/2. + 0.5, H/2. - i + 0.5, -W/(2. * tan(alpha/2.)));
                ray_u.normalize(); 
                Ray ray(Q, ray_u, n1);
                color = color + scene.Color(ray, depth, epsilon);
            }
            color = color/repetition;
            

            //apply some gamma correction
            double gamma_correction = 1./2.2;
            double color_red = pow(color[0], gamma_correction);
            double color_green = pow(color[1], gamma_correction);
            double color_blue = pow(color[2], gamma_correction);


            //create an image  
            image[(i * W + j) * 3 + 0] = std::min(color_red, 255.);
            image[(i * W + j) * 3 + 1] = std::min(color_green, 255.);
            image[(i * W + j) * 3 + 2] = std::min(color_blue, 255.);
        }
    }
    stbi_write_png("image.png", W, H, 3, &image[0], 0);
 
    return 0;
}