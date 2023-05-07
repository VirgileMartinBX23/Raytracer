#ifndef BOUNDINGBOX_H
#define BOUNDINGBOX_H

#include "Vector.h"
#include "Ray.h"



class BoundingBox {
public:
    bool ray_bounding_box_intersection(const Ray& ray) const {
        //according to slides 7 lecture 3
        double tx0 = (B_min[0]-ray.O[0])/ray.u[0];
        double ty0 = (B_min[1]-ray.O[1])/ray.u[1];
        double tz0 = (B_min[2]-ray.O[2])/ray.u[2];

        double tx1 = (B_max[0]-ray.O[0])/ray.u[0];
        double ty1 = (B_max[1]-ray.O[1])/ray.u[1];
        double tz1 = (B_max[2]-ray.O[2])/ray.u[2];

        double tx_min = std::min(tx0, tx1);
        double ty_min = std::min(ty0, ty1);
        double tz_min = std::min(tz0, tz1);

        double tx_max = std::max(tx0, tx1);
        double ty_max = std::max(ty0, ty1);
        double tz_max = std::max(tz0, tz1);
        //according to page 47 of the lecture notes
        double minimum = std::min(tx_max, std::min(ty_max, tz_max));
        double maximum = std::max(tx_min, std::max(ty_min, tz_min));
        if (minimum > maximum){
            return true;
        }
        return false;
    }

    Vector B_min, B_max;
};



#endif