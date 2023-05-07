#ifndef SCENE_H
#define SCENE_H


#include "Vector.h"
#include "Geometry.h"
#include "Ray.h"
#define PI 3.141592653589793238

class Scene {
public:
    void add_sphere(const Sphere* sphere){
        geometries.push_back(sphere);
    }
    void add_mesh(const TriangleMesh* mesh){
        geometries.push_back(mesh);
    }

    bool intersect(const Ray& ray, double &t, Vector &P, Vector &N, int& id, double &n2) {
        t = std::numeric_limits<double>::max();
        bool solution = false;
        for (int i = 0; i < geometries.size(); i++){

            double partial_t, partial_n2;
            Vector partial_P, partial_N;
            
            double partial_solution = (*geometries[i]).Intersection(ray, partial_t, partial_P, partial_N, partial_n2);
            if (partial_solution){
                if (partial_t < t){
                    t = partial_t;
                    P = partial_P;
                    N = partial_N;
                    n2 = partial_n2;
                    id = i;
                    solution = true;

                }
            }
        }
        return solution;
    }


    Vector Color(const Ray& ray, int depth, double epsilon = 0.001){
        Vector color;
        double t, n2; 
        Vector P, N;
        int id;
        double n1 = ray.n1;

        if (depth <= 0){
            return Vector(0,0,0);
        }

        if (intersect(ray, t, P, N, id, n2)){
            if ((*geometries[id]).reflection){
                Vector omega_r = ray.u -2*dot(ray.u, N)*N;
                Ray reflected_ray(P+epsilon*N, omega_r, n1);
                return Color(reflected_ray, depth-1);
            }
            if ((*geometries[id]).refraction) {
				Vector N_refraction = N;
				if (dot(ray.u, N) > 0) { // check if the ray is exiting the sphere
					double n3 = n1;
                    n1 = n2;
                    n2 = n3;
					N_refraction = -1 * N_refraction;
				}
				Vector omega_T = n1/n2 * (ray.u - dot(ray.u, N_refraction)*N_refraction);
				double check = 1 - carre(n1/n2)*(1-carre(dot(ray.u, N_refraction)));
                Vector omega_N = -1 * N_refraction*sqrt(check);
                Vector omega_t = omega_T + omega_N;

				if (check < 0){ //if check < 0 then it's reflective
                    Vector omega_r = ray.u - 2 * dot(ray.u, N) * N;
				    Ray reflected_ray(P - epsilon * N, omega_r, n1);
                    return Color(reflected_ray, depth-1);
                }
				else {
					double k_o = carre(n1-n2)/carre(n1+n2);
					double R = k_o + (1 - k_o)*pow((1 - abs(dot(N_refraction, ray.u))), 5); //reflection coefficient
                    double u = uniform(engine);
					if (u < R){ 
                        Vector omega_r = ray.u - 2 * dot(ray.u, N) * N;
				        Ray reflected_ray(P - epsilon * N, omega_r, n1);
                        return Color(reflected_ray, depth-1);
                    }
					
					Ray refracted_ray(P - epsilon * N_refraction, omega_t, n1);
                    return Color(refracted_ray, depth-1);
				}
            }
            //use to get the color 

            //direct light
            double d = (S - P).norm2(); 
            Vector omega_i = (S - P); 
            omega_i.normalize();

            //the shadow 
            Vector P_prime = P + epsilon*N; //P'
            Ray ray_sha(P_prime, omega_i, n1);
            double t_sha, n2;
            Vector P_sha, N_sha;
            int id_sha;

            if (intersect(ray_sha, t_sha, P_sha, N_sha, id_sha, n2)){
                if (carre(t_sha) < d){ //if the intersection is in the shadow
                    color = Vector(0,0,0);
                }
                else {
                    color = (I*(*geometries[id]).albedo*std::max(0., dot(N, omega_i)))/(4*carre(PI)*d);
                }
            }
            else {
                 color = (I*(*geometries[id]).albedo*std::max(0., dot(N, omega_i)))/(4*carre(PI)*d);

            }

            //indirect light
            
            Vector V = random_cos(N);
            Ray ray_ind(P + epsilon*N, V, n1);
            color = color + (*geometries[id]).albedo*Color(ray_ind, depth-1);

            return color;
        }

        else { //if there is no interaction, there is no color
            color = Vector(0,0,0);
            return color;
        }
    }

    std::vector<const Geometry*> geometries;
    double I;
    Vector S;
    
};
 

#endif