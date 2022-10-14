#pragma once
#include <cassert>
#include <cmath>

namespace krvlib {

bool is_equal(double num_1, double num_2);

struct vector {

    const double x = NAN, y = NAN, z = NAN;


    vector() {} // NAN vector
    
    vector(double x, double y, double z): x(x), y(y), z(z) {}


    bool is_valid() const;

    double mod() const; // Module of vector

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
    
    static vector intersection_point (const line& line_1, const line& line_2);

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
    static vector intersection_point(const segment& segment_1, const segment& segment_2); // if INTERSECTION is POINT return point, else return NAN-vector

    static bool intersection_segment(const segment& segment_1, const segment& segment_2); // If the segments intersect, it returns true, else false

    static INTERSECTION_TYPE intersect(const segment& segment_1, const segment& segment_2);

};



vector intersection_point(const line& _line, const segment& _segment); // if INTERSECTION is POINT return true, else return false



struct plane {

    const vector normal; // Normal vector
    const vector point;  // Radius-vector of a point that lies on the plane


    plane(const vector& normal, const vector& point): normal(normal), point(point) {
        assert(normal != vector(0.0, 0.0, 0.0));
        assert(point.is_valid() && normal.is_valid());
    }


    static line intersection_line(const plane& plane_1, const plane& plane_2); // Liquid function (don't look the definition)

};

bool operator|| (const plane& plane_1, const plane& plane_2); // Check if the planes are parallel

bool operator== (const plane& plane_1, const plane& plane_2); // Check if the planes are equal



struct triangle {

    const vector point_1;
    const vector point_2;
    const vector point_3;


    triangle(const vector& point_1, const vector& point_2, const vector& point_3): point_1(point_1), point_2(point_2), point_3(point_3) {
        assert((point_1 != point_2) || (point_1 != point_3) || (point_2 != point_3));
    }


    plane get_plane() const;

  private:
    segment intersect_on_plane(const line& _line) const; // Works if a triangle and a straight line are in the same plane

  public:    
    static bool intersect(const triangle& triangle_1, const triangle& triangle_2);

};



}