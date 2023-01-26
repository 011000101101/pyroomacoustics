/* 
 * Implementation of the Wall class
 * Copyright (C) 2019  Robin Scheibler, Cyril Cadoux
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * You should have received a copy of the MIT License along with this program. If
 * not, see <https://opensource.org/licenses/MIT>.
 */

#include <iostream>
#include <cmath>

#include "wall.hpp"
#include "geometry.hpp"
#include "common.hpp"

std::shared_ptr<Polygon> Polygon::make_polygon(
        Eigen::Matrix<float, 3, Eigen::Dynamic> corners,
        std::vector<Eigen::Matrix<float, 3, Eigen::Dynamic>> holes
) {
    if(holes.size() > 0) {
        return std::make_shared<PolygonWithHole>(corners, holes);
    }
    else {
        return std::make_shared<SimplePolygon>(corners);
    }
}

std::shared_ptr<Polygon> Polygon::make_polygon(const Polygon &p) {
    if (const PolygonWithHole * derived_ptr = dynamic_cast<const PolygonWithHole *>(&p))
    {
        return std::make_shared<PolygonWithHole>(*derived_ptr);
    }
    else if (const SimplePolygon * derived_ptr = dynamic_cast<const SimplePolygon *>(&p))
    {
        return std::make_shared<SimplePolygon>(*derived_ptr);
    }
}

Polygon::Polygon() {
}

SimplePolygon::SimplePolygon(
        const Eigen::Matrix<float, 3, Eigen::Dynamic> &_corners, const Eigen::Matrix<float, 3, 1> &_origin
) : corners(_corners), origin(_origin) {
    // In 3D things are a little more complicated
    // We need to compute a 2D basis for the plane and find the normal

    // The basis and normal are found by SVD
    Eigen::JacobiSVD<Eigen::Matrix<float,3,Eigen::Dynamic>> svd(corners.colwise() - origin, Eigen::ComputeThinU);

    // The corners matrix should be rank defficient, check the smallest eigen value
    // The rank deficiency is because all the corners are in a 2D subspace of 3D space
    if (svd.singularValues().coeff(2) > libroom_eps)
    {
        throw std::runtime_error("The corners of the polygon do not lie in a plane");
    }

    // The basis is the leading two left singular vectors
    basis.col(0) = svd.matrixU().col(0);
    basis.col(1) = svd.matrixU().col(1);

    // The normal corresponds to the smallest singular value
    normal = svd.matrixU().col(2);

    // Project the 3d corners into 2d plane
    flat_corners = basis.adjoint() * (corners.colwise() - origin);

    // Our convention is that the vertices are arranged counter-clockwise
    // around the normal. In that case, the area computation should be positive.
    // If it is negative, we need to swap the basis.
    float a = area();

    if (a < 0)
    {
        // exchange the other two basis vectors
        basis.rowwise().reverseInPlace();
        flat_corners.colwise().reverseInPlace();
    }

    // Now the normal is computed as the cross product of the two basis vectors
    normal = cross(basis.col(0), basis.col(1));
}

PolygonWithHole::PolygonWithHole(
        const Eigen::Matrix<float, 3, Eigen::Dynamic> &_corners,
        const std::vector<Eigen::Matrix<float, 3, Eigen::Dynamic>> &_holes
) : outer_polygon(SimplePolygon(_corners)) {
    for (unsigned i=0; i < _holes.size(); i++) {
        inner_polygons.push_back(SimplePolygon(_holes[i]));
    }
}

PolygonWithHole::PolygonWithHole(const PolygonWithHole &p) : outer_polygon(SimplePolygon(p.get_outer_polygon())) {
    std::vector<SimplePolygon> original_inner_polygons = p.get_inner_polygons();
    for (unsigned i=0; i < original_inner_polygons.size(); i++) {
        inner_polygons.push_back(SimplePolygon(original_inner_polygons[i]));
    }
}

float Wall2D::area() const
{
  return (corners.col(1) - corners.col(0)).norm();
}

float SimplePolygon::area() const
{
    return area_2d_polygon(flat_corners);
}

float PolygonWithHole::area() const
{
    float inner_areas = 0.;
    for (unsigned i=0; i < inner_polygons.size(); i++) {
        inner_areas += inner_polygons[i].area();
    }
    return outer_polygon.area() - inner_areas;
}

int SimplePolygon::intersection(
        const Eigen::Matrix<float,3,1> &p1,
        const Eigen::Matrix<float,3,1> &p2,
        Eigen::Ref<Eigen::Matrix<float,3,1>> intersection
) const
{
    /*
      Computes the intersection between a line segment and a polygon surface in 3D.
      This function computes the intersection between a line segment (defined
      by the coordinates of two points) and a surface (defined by an array of
      coordinates of corners of the polygon and a normal vector)
      If there is no intersection, None is returned.
      If the segment belongs to the surface, None is returned.
      Two booleans are also returned to indicate if the intersection
      happened at extremities of the segment or at a border of the polygon,
      which can be useful for limit cases computations.

      a1: (array size 3) coordinates of the first endpoint of the segment
      a2: (array size 3) coordinates of the second endpoint of the segment
      corners: (array size 3xN, N>2) coordinates of the corners of the polygon
      normal: (array size 3) normal vector of the surface
      intersection: (array size 3) store the intersection point

      :returns:
             -1 if there is no intersection
              0 if the intersection striclty between the segment endpoints and in the polygon interior
              1 if the intersection is at endpoint of segment
              2 if the intersection is at boundary of polygon
              3 if both the above are true
      */

    int ret1, ret2, ret = 0;

    ret1 = intersection_3d_segment_plane(p1, p2, origin, normal, intersection);

    if (ret1 == -1)
        return -1;  // there is no intersection

    if (ret1 == 1)  // intersection at endpoint of segment
        ret = 1;

    /* project intersection into plane basis */
    Eigen::Vector2f flat_intersection = basis.adjoint() * (intersection - origin);

    /* check in flatland if intersection is in the polygon */
    ret2 = is_inside_2d_polygon(flat_intersection, flat_corners);

    if (ret2 < 0)  // intersection is outside of the wall
        return -1;

    if (ret2 == 1) // intersection is on the boundary of the wall
        ret |= 2;

    return ret;  // no intersection
}

int PolygonWithHole::intersection(
        const Eigen::Matrix<float,3,1> &p1,
        const Eigen::Matrix<float,3,1> &p2,
        Eigen::Ref<Eigen::Matrix<float,3,1>> intersection
) const
{
    int ret1 = outer_polygon.intersection(p1, p2, intersection);
    // if there is no intersection with the outer polygon, there can be no intersection overall
    if (ret1 == -1) {
        return ret1;
    }
    // if the intersection lies on the boundary of the outer polygon, it is not affected by inner polygons, as we do
    // not allow holes to coincide with the outer boundary
    else if (ret1 >= 2) {
        return ret1;
    }
    else {
        // check intersection with each of the holes
        for (unsigned i=0; i < inner_polygons.size(); i++) {
            int ret2 = inner_polygons[i].intersection(p1, p2, intersection);
            // if the intersection falls within one hole, it does not intersect the total polygon and we can return
            if (ret2 == 0 || ret2 == 1) {
                return -1;
            }
            // the boundary of holes is treated like the boundary of the total polygon
            else if (ret2 >= 2) {
                return ret2;
            }
            // else continue with the next hole
        }
        // if no intersection with any hole, return the original intersection value with the outer polygon
        return ret1;
    }
}

bool SimplePolygon::same_as(const std::shared_ptr<Polygon> that) const
{
    /*
    Checks if two walls are the same, based on their corners of the walls.
    Be careful : it will return true for two identical walls that belongs
    to two different rooms !
    */
    if (const std::shared_ptr<const SimplePolygon> other = std::dynamic_pointer_cast<const SimplePolygon>(that))
    {
        // Not the same number of corners
        if (corners.cols() != other->get_corners().cols())
        {
            return false;
        }
        // check corner coordinates
        return (corners - other->get_corners()).cwiseAbs().sum() == 0.;
    }
    else
    {
        // The two walls are not of the same geometry: One has holes, the other does not
        return false;
    }
}

bool PolygonWithHole::same_as(const std::shared_ptr<Polygon> that) const
{
    /*
    Checks if two walls are the same, based on their corners of the walls.
    Be careful : it will return true for two identical walls that belongs
    to two different rooms !
    */
    if (const std::shared_ptr<const PolygonWithHole> other = std::dynamic_pointer_cast<const PolygonWithHole>(that))
    {
        if (!outer_polygon.same_as(std::make_shared<SimplePolygon>(other->get_outer_polygon())))
        {
            // outer does not match
            return false;
        }
        if (inner_polygons.size() != other->get_inner_polygons().size())
        {
            // number of holes does not match
            return false;
        }
        // both vectors have the same size and can thus be indexed safely by the same variable
        for (unsigned int i = 0; i < inner_polygons.size(); i++)
        {
            // TODO assumes holes to be sorted
            if (!inner_polygons[i].same_as(std::make_shared<SimplePolygon>(other->get_inner_polygons()[i])))
            {
                // a hole does not match
                return false;
            }
        }
        // everything matches
        return true;
    }
    else
    {
        // The two walls are not of the same geometry: One has holes, the other does not
        return false;
    }
}

float Wall3D::area() const
{
    return wall_geometry->area();
}

template<size_t D>
void Wall<D>::init()
{
  // compute transmission coefficients from absorption
  energy_reflection.resize(absorption.size());
  energy_reflection = 1.f - absorption;
  transmission.resize(absorption.size());
  transmission = energy_reflection.sqrt();

  if (absorption.size() != scatter.size())
  {
    throw std::runtime_error("The number of absorption and scattering coefficients is different");
  }
}

template<size_t D>
Wall<D>::Wall(
    const Eigen::Matrix<float,D,Eigen::Dynamic> &_corners,
    const Eigen::ArrayXf &_absorption,
    const Eigen::ArrayXf &_scatter,
    const std::string &_name
    ) : corners(_corners), absorption(_absorption), scatter(_scatter), name(_name)
{
    init();
}

//template<>
Wall2D::Wall2D(
    const Eigen::Matrix<float,2,Eigen::Dynamic> &_corners,
    const Eigen::ArrayXf &_absorption,
    const Eigen::ArrayXf &_scatter,
    const std::string &_name
    )
  : Wall<2>(_corners, _absorption, _scatter, _name)
{

  // Pick one of the corners as the origin of the wall
  origin = corners.col(0);

  // compute normal (difference of 2 corners, swap x-y, change 1 sign)
  normal.coeffRef(0) = corners.coeff(1,1) - corners.coeff(1,0);
  normal.coeffRef(1) = corners.coeff(0,0) - corners.coeff(0,1);
  normal = normal.normalized();
}

//template<>
Wall3D::Wall3D(
    const Eigen::Matrix<float,3,Eigen::Dynamic> &_corners,
    const std::vector<Eigen::Matrix<float, 3, Eigen::Dynamic>> &_holes,
    const Eigen::ArrayXf &_absorption,
    const Eigen::ArrayXf &_scatter,
    const std::string &_name
    )
  : Wall<3>(_corners, _absorption, _scatter, _name), wall_geometry(Polygon::make_polygon(_corners, _holes))
{

    // shadow members of wall_geometry to allow direct access from outside for backward compatibility
    // shadowed members are constant, therefore these copies are always valid

    // Pick the origin as the first corner of the outer polygon
    origin = wall_geometry->get_origin();
//    corners = wall_geometry->get_corners();  // already set by Wall<3> constructor
    //basis = wall_geometry->get_basis();
//    normal = wall_geometry->get_normal();
    //flat_corners = wall_geometry->get_flat_corners();
    normal = wall_geometry->get_normal();
}

int Wall2D::intersection(
    const Eigen::Matrix<float,2,1> &p1,
    const Eigen::Matrix<float,2,1> &p2,
    Eigen::Ref<Eigen::Matrix<float,2,1>> intersection
    ) const
{
  return intersection_2d_segments(p1, p2, corners.col(0), corners.col(1), intersection);
}

//template<>
int Wall3D::intersection(
    const Eigen::Matrix<float,3,1> &p1,
    const Eigen::Matrix<float,3,1> &p2,
    Eigen::Ref<Eigen::Matrix<float,3,1>> intersection
    ) const
{
  /*
    Computes the intersection between a line segment and a polygon surface in 3D.
    This function computes the intersection between a line segment (defined
    by the coordinates of two points) and a surface (defined by an array of
    coordinates of corners of the polygon and a normal vector)
    If there is no intersection, None is returned.
    If the segment belongs to the surface, None is returned.
    Two booleans are also returned to indicate if the intersection
    happened at extremities of the segment or at a border of the polygon,
    which can be useful for limit cases computations.

    a1: (array size 3) coordinates of the first endpoint of the segment
    a2: (array size 3) coordinates of the second endpoint of the segment
    corners: (array size 3xN, N>2) coordinates of the corners of the polygon
    normal: (array size 3) normal vector of the surface
    intersection: (array size 3) store the intersection point

    :returns:
           -1 if there is no intersection
            0 if the intersection striclty between the segment endpoints and in the polygon interior
            1 if the intersection is at endpoint of segment
            2 if the intersection is at boundary of polygon
            3 if both the above are true
    */

  return wall_geometry->intersection(p1, p2, intersection);
}

template<size_t D>
int Wall<D>::intersects(const Vectorf<D> &p1, const Vectorf<D> &p2) const
{
  Vectorf<D> v;
  return intersection(p1, p2, v);
}

template<size_t D>
int Wall<D>::reflect(const Vectorf<D> &p, Eigen::Ref<Vectorf<D>> p_reflected) const
{
  /*
   * Reflects point p across the wall 
   *
   * wall: a wall object (2d or 3d)
   * p: a point in space
   * p_reflected: a pointer to a buffer large enough to receive
   *              the location of the reflected point
   *
   * Returns: 1 if reflection is in the same direction as the normal
   *          0 if the point is within tolerance of the wall
   *         -1 if the reflection is in the opposite direction of the normal
   */

  // projection onto normal axis
  float distance_wall2p = normal.adjoint() * (origin - p);

  // compute reflected point
  p_reflected = p + 2 * distance_wall2p * normal;

  if (distance_wall2p > libroom_eps)
    return 1;
  else if (distance_wall2p < -libroom_eps)
    return -1;
  else
    return 0;
}


/* checks on which side of a wall a point is */
template<size_t D>
int Wall<D>::side(const Vectorf<D> &p) const
{
  // Essentially, returns the sign of the inner product with the normal vector
  float ip = (p - origin).adjoint() * normal;

  if (ip > libroom_eps)
    return 1;
  else if (ip < -libroom_eps)
    return -1;
  else
    return 0;
}

bool Wall2D::same_as(const Wall & that) const
{
  /*
  Checks if two walls are the same, based on their corners of the walls.
  Be careful : it will return true for two identical walls that belongs
  to two different rooms !
  */

  if (that.dim != 2)
  {
    std::cerr << "The two walls are not of the same dimensions !" << std::endl;
    return false;
  }

  // Not the same number of corners
  if (corners.cols() != that.corners.cols())
  {
    return false;
  }

  return (corners - that.corners).cwiseAbs().sum() == 0.;
}

bool Wall3D::same_as(const Wall & that) const
{
    /*
    Checks if two walls are the same, based on their corners of the walls.
    Be careful : it will return true for two identical walls that belongs
    to two different rooms !
    */

    if (that.dim != 3)
    {
        std::cerr << "The two walls are not of the same dimensions !" << std::endl;
        return false;
    }
    // downcast to Wall3D to be able to access wall_geometry
    // TODO use static_cast?
    const Wall3D * that_specific = dynamic_cast<const Wall3D *>(&that);
    // check if the geometry is the same
    return wall_geometry->same_as(that_specific->wall_geometry);
}


Eigen::Matrix<float, 2, Eigen::Dynamic> Wall2D::get_corners() const
{
    return corners;
}

Eigen::Matrix<float, 3, Eigen::Dynamic> Wall3D::get_corners() const
{
    return  wall_geometry->get_corners();
}


std::vector<Eigen::Matrix<float, 2, Eigen::Dynamic>> Wall2D::get_holes() const
{
    return std::vector<Eigen::Matrix<float, 2, Eigen::Dynamic>>();
}

std::vector<Eigen::Matrix<float, 3, Eigen::Dynamic>> Wall3D::get_holes() const
{
    return  wall_geometry->get_holes();
}

template<size_t D>
Vectorf<D> Wall<D>::normal_reflect(
    const Vectorf<D> &start,
    const Vectorf<D> &hit_point,
    float length) const
{
	  
  /* This method computes the reflection of one point with respect to
   a precise hit_point on a wall. Also, the distance between the
   wall hit point and the reflected point is defined by the 'length'
   parameter.
   This method computes the reflection of point 'start' across the normal
   to the wall through 'hit_point'.
    
   start: (array size 2 or 3) defines the point to be reflected
   hit_point: (array size 2 or 3) defines a point on a wall that will
     serve as the reference point for the reflection
   wall_normal: (array size 2 or 3) defines the normal of the reflecting
     wall. It will be used as if it was anchored at hit_point
   length : the desired distance between hit_point and the reflected point
   
   :returns: an array of size 2 or 3 representing the reflected point
   */

  Vectorf<D> incident = (hit_point - start).normalized();
  // return hit_point + length * normal_reflect(incident);
  return hit_point + length * (incident - normal * 2 * incident.dot(normal));
}

template<size_t D>
Vectorf<D> Wall<D>::normal_reflect(const Vectorf<D> &incident) const
{
  /*
   * Same as the previous function, but works on a direction vector instead
   */
  return incident - normal * 2 * incident.dot(normal);
}

template<size_t D>
float Wall<D>::cosine_angle(
    const Vectorf<D> &p) const
{
    /*
    Compute the cosine angle between the surface normal and a given vector.
    */

    return p.dot(normal) / p.norm();

}


