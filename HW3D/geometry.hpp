#pragma once
#include <cassert>
#include <cmath>
#include <utility>
#include <algorithm>

namespace krvlib {

bool is_equal(double num_1, double num_2);


struct vector {

    const double x = NAN, y = NAN, z = NAN;


    vector() {} // NAN vector
    
    vector(double x, double y, double z): x(x), y(y), z(z) {}


    bool is_valid() const;

    double mod() const; // Module of vector

    vector normal() const; //return normal vector

    static vector vect_mp(const vector& vector_1, const vector& vector_2); // Returns a vector parallel to the vector multiplication


    vector operator- () const; // Don't have protection from notfinity
};

bool operator== (const vector& vector_1, const vector& vector_2); // Don't have protection from notfinity

bool operator!= (const vector& vector_1, const vector& vector_2); // Don't have protection from notfinity

vector operator- (const vector& vector_1, const vector& vector_2); // Don't have protection from notfinity

vector operator+ (const vector& vector_1, const vector& vector_2); // Don't have protection from notfinity

vector operator* (const vector& vect, const double coef); // Don't have protection from notfinity

double operator* (const vector& vector_1, const vector& vector_2); //Scalar multiplication

bool operator|| (const vector& vector_1, const vector& vector_2);  // Checks if the vectors are parallel



struct line {
    
    const vector vec_d = vector(); // Guiding vector
    const vector point = vector(); // Radius-vector of a point through which the line passes

    
    line() {}

    line(const vector& vec_d, const vector& point): vec_d(vec_d), point(point) {\
        assert(point.is_valid() && vec_d.is_valid());
    }


    bool is_valid() const;

    bool belong(const vector& point) const;

    static bool intersect_once_on_plane(const line& line_1, const line& line_2);

  private:
    static vector intersection_point_function(const line& line_1, const line& line_2); // Liquid function (don't look the definition)

  public:
    
    static vector intersection_point (const line& line_1, const line& line_2);//Works if lines are in the same plane

};

bool operator|| (const line& line_1, const line& line_2); // Checks if the lines are parallel

bool operator== (const line& line_1, const line& line_2); // Checks if the lines are equal



struct segment {

    const vector point_1 = vector(); // Beginning of the segment
    const vector point_2 = vector(); // End of the segment

    enum class INTERSECTION_TYPE {
        NONE,
        POINT,
        SEGMENT,
    };



    segment() {}

    segment(const vector& point_1, const vector& point_2): point_1(point_1), point_2(point_2) {
        assert(point_1.is_valid() && point_2.is_valid());
    }


    bool is_valid() const;

    bool is_point() const;

    line get_line() const;

    bool belong(const vector& point) const;

  private:
    static bool boundary_points_equal(const segment& segment_1, const segment& segment_2);

  public:
    static vector intersection_point(const segment& segment_1, const segment& segment_2); // if INTERSECTION is POINT return point, else return NAN-vector (same plane)

    static bool intersection_segment(const segment& segment_1, const segment& segment_2); // If the segments intersect, it returns true, else false (same plane)

    static INTERSECTION_TYPE intersect(const segment& segment_1, const segment& segment_2);

};



vector intersection_point(const line& _line, const segment& _segment); // if INTERSECTION is POINT return true, else return false (same plane)



struct plane {

    const vector normal; // Normal vector
    const vector point;  // Radius-vector of a point that lies on the plane


    plane(const vector& normal, const vector& point): normal(normal), point(point) {
        assert(point.is_valid() && normal.is_valid());
    }


    bool is_valid() const;

    static line intersection_line(const plane& plane_1, const plane& plane_2); // Liquid function (don't look the definition)

};

bool operator|| (const plane& plane_1, const plane& plane_2); // Check if the planes are parallel

bool operator== (const plane& plane_1, const plane& plane_2); // Check if the planes are equal


struct AABB_box {
    const double x_min, x_max;
    const double y_min, y_max;
    const double z_min, z_max;


    AABB_box(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max): x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max), z_min(z_min), z_max(z_max) {}
    // AABB_box(const triangle& _triangle): x_min(std::min({_triangle.point_1.x, _triangle.point_2.x, _triangle.point_3.x})), x_max(std::max({_triangle.point_1.x, _triangle.point_2.x, _triangle.point_3.x})),
    // y_min(std::min({_triangle.point_1.y, _triangle.point_2.y, _triangle.point_3.y})), y_max(std::max({_triangle.point_1.y, _triangle.point_2.y, _triangle.point_3.y})),
    // z_min(std::min({_triangle.point_1.z, _triangle.point_2.z, _triangle.point_3.z})), z_max(std::max({_triangle.point_1.z, _triangle.point_2.z, _triangle.point_3.z}))
    // {}

    static bool intersect(const AABB_box& box_1, const AABB_box& box_2);
};


struct triangle {

    const vector point_1;
    const vector point_2;
    const vector point_3;


    triangle(const vector& point_1, const vector& point_2, const vector& point_3): point_1(point_1), point_2(point_2), point_3(point_3) {
    }


    plane get_plane() const;

    AABB_box get_AABB_box() const;

  private:
    segment intersect_on_plane(const line& _line) const; // Works if a triangle and a straight line are in the same plane

    plane plane_getter_for_special_case_intersects() const;

    static std::pair<plane, plane> plane_getter_for_intersect(const triangle& triangle_1, const triangle& triangle_2);

    static bool intersect_accurate_extra(const triangle& triangle_1, const triangle& triangle_2);

  public:    
    static bool intersect(const triangle& triangle_1, const triangle& triangle_2);

};



}