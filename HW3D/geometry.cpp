#include <cassert>
#include <limits>
#include "geometry.hpp"

namespace krvlib {

////////////////////////////////////////////////////vectors///////////////////////////////////////////////////////
bool is_equal(double num_1, double num_2) {
    if (!(std::isfinite(num_1)) || !(std::isfinite(num_2))) {
        return false;
    }
    double plus = std::abs(num_1 + num_2), minus = std::abs(num_1 - num_2);
    if (minus < (std::numeric_limits<double>::epsilon())*plus || minus < static_cast<double>(std::numeric_limits<double>::min())) {
        return true;
    }
    return false;
}

bool vector::is_valid() const {
    return (std::isfinite(x) && std::isfinite(y) && std::isfinite(z));
}

double vector::mod() const {
    assert(is_valid());
    return sqrt(x*x + y*y + z*z);
}

vector vector::vect_mp(const vector& vector_1, const vector& vector_2) {
    assert(vector_1.is_valid() && vector_2.is_valid());
    return vector(vector_1.y*vector_2.z - vector_1.z*vector_2.y, vector_1.z*vector_2.x - vector_1.x*vector_2.z, vector_1.x*vector_2.y - vector_1.y*vector_2.x);
}

vector vector::operator- () const {
    return vector(-x, -y, -z);
}

bool operator== (const vector& vector_1, const vector& vector_2) {
    return (is_equal(vector_1.x, vector_2.x) && is_equal(vector_1.y, vector_2.y) && is_equal(vector_1.z, vector_2.z));
}

bool operator!= (const vector& vector_1, const vector& vector_2) {
    return !(vector_1 == vector_2);
}

vector operator+ (const vector& vector_1, const vector& vector_2) {
    return vector(vector_1.x + vector_2.x, vector_1.y + vector_2.y, vector_1.z + vector_2.z);
}

vector operator- (const vector& vector_1, const vector& vector_2) {
    return vector(vector_1 + (-vector_2));
}

vector operator* (const vector& vect, const double coef) {
    return vector(vect.x*coef, vect.y*coef, vect.z*coef);
}

double operator* (const vector& vector_1, const vector& vector_2) {
    assert(vector_1.is_valid() && vector_2.is_valid());
    return vector_1.x*vector_2.x + vector_1.y*vector_2.y + vector_1.z*vector_2.z;
}

bool operator|| (const vector& vector_1, const vector& vector_2) {
    return is_equal(std::abs(vector_1*vector_2), vector_1.mod()*vector_2.mod());
}

///////////////////////////////////////////////////////////lines//////////////////////////////////////////////////////////

bool line::is_valid() const {
    return ((vec_d != vector(0.0, 0.0, 0.0)) && (vec_d.is_valid()) && (point.is_valid()));
}

bool line::belong(const vector& point) const {
    assert(is_valid() && point.is_valid());
    return (vec_d||(this->point - point));
}

bool line::intersect_once_on_plane(const line& line_1, const line& line_2) {//PROoooooooooooooooooBLEMMMMMMMMMM
    return !( line_1||line_2 );
}

vector line::intersection_point_function(const line& line_1, const line& line_2) {
    // Liquid moment
    // A system of equations of the form: v*t + x = u*s + y
    //          
    //          v1*t + x1 = u1*s + y1
    //          v2*t + x2 = u2*s + y2
    //          v3*t + x3 = u3*s + y3
    //
    //          1, 2 and 3 are the x, y and z coordinates
    //
    // It is guaranteed that the system of three equations has solutions t and s
    assert(line_1.is_valid() && line_2.is_valid());
    double v[3] = {line_1.vec_d.x, line_1.vec_d.y, line_1.vec_d.z};
    double x[3] = {line_1.point.x, line_1.point.y, line_1.point.z};
    double u[3] = {line_2.vec_d.x, line_2.vec_d.y, line_2.vec_d.z};
    double y[3] = {line_2.point.x, line_2.point.y, line_2.point.z};

    if (!is_equal(v[0], 0.0)) { // Case x
        double v1 = v[0], x1 = x[0], u1 = u[0], y1 = y[0];
        if (!is_equal(v1*u[1] - v[1]*u1, 0.0)) {
            double v2 = v[1], x2 = x[1], u2 = u[1], y2 = y[1];
            double s = (v1*(x2 - y2) - v2*(x1 - y1))/(v1*u2 - v2*u1);
            return (line_2.vec_d)*s + line_2.point;
        }
        if (!is_equal(v1*u[2] - v[2]*u1, 0.0)) {
            double v2 = v[2], x2 = x[2], u2 = u[2], y2 = y[2];
            double s = (v1*(x2 - y2) - v2*(x1 - y1))/(v1*u2 - v2*u1);
            return (line_2.vec_d)*s + line_2.point;
        }
    }
    else if (!is_equal(v[1], 0.0)) { // Case y
        double v1 = v[1], x1 = x[1], u1 = u[1], y1 = y[1];
        if (!is_equal(v1*u[0] - v[0]*u1, 0.0)) {
            double v2 = v[0], x2 = x[0], u2 = u[0], y2 = y[0];
            double s = (v1*(x2 - y2) - v2*(x1 - y1))/(v1*u2 - v2*u1);
            return (line_2.vec_d)*s + line_2.point;
        }
        if (!is_equal(v1*u[2] - v[2]*u1, 0.0)) {
            double v2 = v[2], x2 = x[2], u2 = u[2], y2 = y[2];
            double s = (v1*(x2 - y2) - v2*(x1 - y1))/(v1*u2 - v2*u1);
            return (line_2.vec_d)*s + line_2.point;
        }
    }
    else if (!is_equal(v[2], 0.0)) { // Case z
        double v1 = v[2], x1 = x[2], u1 = u[2], y1 = y[2];
        if (!is_equal(v1*u[0] - v[0]*u1, 0.0)) {
            double v2 = v[0], x2 = x[0], u2 = u[0], y2 = y[0];
            double s = (v1*(x2 - y2) - v2*(x1 - y1))/(v1*u2 - v2*u1);
            return (line_2.vec_d)*s + line_2.point;
        }
        if (!is_equal(v1*u[1] - v[1]*u1, 0.0)) {
            double v2 = v[1], x2 = x[1], u2 = u[1], y2 = y[1];
            double s = (v1*(x2 - y2) - v2*(x1 - y1))/(v1*u2 - v2*u1);
            return (line_2.vec_d)*s + line_2.point;
        }
    }
    assert(false);
}

vector line::intersection_point (const line& line_1, const line& line_2) {
    assert(line_1.is_valid() && line_2.is_valid());
    if (!line::intersect_once_on_plane(line_1, line_2)) {
        return vector();
    }
    if (line_1.point == line_2.point) {
        return line_1.point;
    }

    return intersection_point_function(line_1, line_2);
}

bool operator|| (const line& line_1, const line& line_2) {
    return line_1.vec_d||line_2.vec_d;
}

bool operator== (const line& line_1, const line& line_2) {
    return ((line_1.vec_d||line_2.vec_d) && (line_1.vec_d||(line_1.point - line_2.point)));
}

///////////////////////////////////////////////////////segments////////////////////////////////////////////////////////////

bool segment::is_valid() const {
    return ((point_1.is_valid() && point_2.is_valid()));  
}

bool segment::is_point() const {
    return (point_1 == point_2);
}

line segment::get_line() const {
    return line(point_2 - point_1, point_1);
}

bool segment::belong(const vector& point) const {
    if (!(point.is_valid())) {
        return false;
    }
    return ((is_equal((point_1 - point).mod() + (point_2 - point).mod(), (point_2 - point_1).mod())));
}

bool segment::boundary_points_equal(const segment& segment_1, const segment& segment_2) {
    const vector& p11 = segment_1.point_1, p12 = segment_1.point_2, p21 = segment_2.point_1, p22 = segment_2.point_2;
    return ((p11 == p21) || (p11 == p22) || (p12 == p21) || (p12 == p22));
}

vector segment::intersection_point (const segment& segment_1, const segment& segment_2) {
    if (segment_1.is_point()) {
        if (segment_2.belong(segment_1.point_1)) {
            return segment_1.point_1;
        }
        return vector();
    }
    if (segment_2.is_point()) {
        if (segment_1.belong(segment_2.point_1)) {
            return segment_2.point_1;
        }
        return vector();
    }
    const line& line_1 = segment_1.get_line(), line_2 = segment_2.get_line();
    if (line_1 == line_2) {
        const vector& s1p1 = segment_1.point_1;
        const vector& s1p2 = segment_1.point_2;
        const vector& s2p1 = segment_2.point_1;
        const vector& s2p2 = segment_2.point_2;
        if ((s1p1 == s2p1) && (((s1p2 - s1p1)*(s2p2 - s2p1)) < 0.0)) {
            return s1p1;
        }
        if ((s1p2 == s2p1) && (((s1p1 - s1p2)*(s2p2 - s2p1)) < 0.0)) {
            return s1p2;
        }
        if ((s1p1 == s2p2) && (((s1p2 - s1p1)*(s2p1 - s2p2)) < 0.0)) {
            return s1p1;
        }
        if ((s1p2 == s2p2) && (((s1p1 - s1p2)*(s2p1 - s2p2)) < 0.0)) {
            return s1p2;
        }
    }
    vector point = line::intersection_point(line_1, line_2);
    if (segment_1.belong(point) && segment_2.belong(point)) {
        return point;
    }
    return vector();
}

bool segment::intersection_segment(const segment& segment_1, const segment& segment_2) {
    if (!(segment_1.is_valid() && segment_2.is_valid())) {
        return false;
    }
    const vector& s1p1 = segment_1.point_1;
    const vector& s1p2 = segment_1.point_2;
    const vector& s2p1 = segment_2.point_1;
    const vector& s2p2 = segment_2.point_2;
    if (segment_2.belong(s1p1) || segment_2.belong(s1p2) || segment_1.belong(s2p1) || segment_1.belong(s2p2)) {
        return true;
    }
    return false;
}

segment::INTERSECTION_TYPE segment::intersect(const segment& segment_1, const segment& segment_2) {
    assert(segment_1.is_valid() && segment_2.is_valid());
    if (intersection_point(segment_1, segment_2).is_valid()) {
        return INTERSECTION_TYPE::POINT;
    }
    if (intersection_segment(segment_1, segment_2)) {
        return INTERSECTION_TYPE::SEGMENT;
    }
    return INTERSECTION_TYPE::NONE;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector intersection_point(const line& _line, const segment& _segment) {
    assert(_line.is_valid());
    if (!(_segment.is_valid())) {
        return vector();
    }
    vector point = line::intersection_point(_line, _segment.get_line());
    if (_segment.belong(point)) {
        return point;
    }
    return vector();
}

////////////////////////////////////////////////////////planes/////////////////////////////////////////////////////////////

line plane::intersection_line(const plane& plane_1, const plane& plane_2) {
    if ( plane_1||plane_2 ) {
        return line();
    }
    // Liquid moment ///////////////////////////////////////////////
    // A system of equations of the form: a1*x + b1*y + c1*z + d1 = 0
    //          
    //          a1*x + b1*y + c1*z = d1 // The first plane
    //          a2*x + b2*y + c2*z = d2 // The second plane
    //
    ////////////////////////////////////////////////////////////////
    vector vec_d = vector::vect_mp(plane_1.normal, plane_2.normal);

    if (plane_1.point == plane_2.point) {
        return line(vec_d, plane_1.point);
    }

    double a1 = plane_1.normal.x, b1 = plane_1.normal.y, c1 = plane_1.normal.z, d1 = (plane_1.normal*plane_1.point);
    double a2 = plane_2.normal.x, b2 = plane_2.normal.y, c2 = plane_2.normal.z, d2 = (plane_2.normal*plane_2.point);
    
    double delta = a1*b2 - b1*a2; // Cramer's rule
    if (!is_equal(delta, 0.0)) { // z = 0
        double delta_x = d1*b2 - b1*d2, delta_y = a1*d2 - d1*a2;
        return line(vec_d, vector(delta_x/delta, delta_y/delta, 0.0));
    }

    delta = b1*c2 - c1*b2;
    if (!is_equal(delta, 0.0)) { // x = 0
        double delta_y = d1*c2 - c1*d2, delta_z = b1*d2 - d1*b2;
        return line(vec_d, vector(0.0, delta_y/delta, delta_z/delta));
    }

    delta = a1*c2 - c1*a2;
    if (!is_equal(delta, 0.0)) { // y = 0
        double delta_x = d1*c2 - c1*d2, delta_z = a1*d2 - d1*a2;
        return line(vec_d, vector(delta_x/delta, 0.0, delta_z/delta));
    }
    assert(false);
}

bool operator|| (const plane& plane_1, const plane& plane_2) {
    return plane_1.normal||plane_2.normal;
}

bool operator== (const plane& plane_1, const plane& plane_2) {
    if ( plane_1.normal||plane_2.normal ) {
        if (is_equal(((plane_1.point - plane_2.point)*plane_1.normal), 0.0)) {
            return true;
        }
    }
    return false;
}

//////////////////////////////////////////////triangles////////////////////////////////////////////////////

plane triangle::get_plane() const {
    return plane(vector::vect_mp(point_3 - point_1, point_2 - point_1), point_1);
}

segment triangle::intersect_on_plane(const line& _line) const {
    assert(_line.is_valid());
    segment sgm_1(point_1, point_2), sgm_2(point_2, point_3), sgm_3(point_3, point_1);
    vector p_1 = intersection_point(_line, sgm_1);
    vector p_2 = intersection_point(_line, sgm_2);
    if ((!p_1.is_valid()) && (!p_2.is_valid())) {
        return segment();
    }
    vector p_3 = intersection_point(_line, sgm_3);

    if (p_1.is_valid() && p_2.is_valid() && p_3.is_valid()) {
        if (p_1 == p_2) {
            return segment(p_2, p_3);
        }
        if (p_2 == p_3) {
            return segment(p_1, p_3);
        }
        if (p_3 == p_1) {
            return segment(p_1, p_2);
        }
    }
    if (!p_1.is_valid()) {
        return segment(p_2, p_3);
    }
    if (!p_2.is_valid()) {
        return segment(p_1, p_3);
    }
    if (!p_3.is_valid()) {
        return segment(p_1, p_2);
    }
    assert(false);
}

bool triangle::intersect(const triangle& triangle_1, const triangle& triangle_2) {
    const plane& plane_1 = triangle_1.get_plane();
    const plane& plane_2 = triangle_2.get_plane();
    if (plane_1 == plane_2) {
        const vector& point1_1 = triangle_1.point_1, point1_2 = triangle_1.point_2, point1_3 = triangle_1.point_3;
        const vector& point2_1 = triangle_2.point_1, point2_2 = triangle_2.point_2, point2_3 = triangle_2.point_3;
        segment sgm1[3] = {segment(point1_2, point1_1), segment(point1_3, point1_2), segment(point1_1, point1_3)};
        segment sgm2[3] = {segment(point2_2, point2_1), segment(point2_3, point2_2), segment(point2_1, point2_3)};
        for (int i = 0; i < 3; i++) {
            if (segment::intersection_segment(sgm1[i], triangle_2.intersect_on_plane((sgm1[i]).get_line()))) {
                return true;
            }
        }
        for (int i = 0; i < 3; i++) {
            if (segment::intersection_segment(sgm2[i], triangle_1.intersect_on_plane((sgm2[i]).get_line()))) {
                return true;
            }
        }
        return false;
    }
    if (plane_1 || plane_2) {
        return false;
    }

    line intersection_line = plane::intersection_line(plane_1, plane_2);

    return segment::intersection_segment(triangle_1.intersect_on_plane(intersection_line), triangle_2.intersect_on_plane(intersection_line));
}



}