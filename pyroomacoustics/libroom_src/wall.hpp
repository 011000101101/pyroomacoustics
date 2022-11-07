/* 
 * Definition of the Wall class used in libroom core of pyroomacoustics
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
#ifndef __CWALL_H__
#define __CWALL_H__

#include <string>
#include <Eigen/Dense>

extern float libroom_eps;

#define WALL_ISECT_NONE        -1  // if there is no intersection
#define WALL_ISECT_VALID        0  // if the intersection striclty between the segment endpoints and in the polygon interior
#define WALL_ISECT_VALID_ENDPT  1  // if the intersection is at endpoint of segment
#define WALL_ISECT_VALID_BNDRY  2  // if the intersection is at boundary of polygon
#define ENDPOINT_BOUNDARY       3  // if both the above are true

//enum class Isect {  // The different cases for intersections
//    NONE = -1,  // - There is no intersection
//    VALID = 0,  // - There is a valid intersection
//    ENDPT = 1,  // - The intersection is on the endpoint of the segment
//    BNDRY = 2   // - The intersection is on the boundary of the wall
//};

class Polygon
{
    private:

    public:

        //factory method
        static Polygon *make_polygon(
                Eigen::Matrix<float, 3, Eigen::Dynamic> corners,
                std::vector<Eigen::Matrix<float, 3, Eigen::Dynamic>> holes
        );

        // Constructor
        Polygon();

        virtual float area() const = 0;  // compute the area of the wall
        virtual Eigen::Matrix<float, 3, 1> get_origin() const = 0;
        virtual int intersection(  // compute the intersection of line segment (p1 <-> p2) with wall
                const Vectorf<3> &p1,
                const Vectorf<3> &p2,
                Eigen::Ref<Vectorf<3>> intersection
        ) const = 0;

//        virtual int reflect(
//                const Vectorf<3> &p,
//                Eigen::Ref<Vectorf<3>> p_reflected
//        ) const;
//        virtual int side(const Vectorf<3> &p) const;
//        virtual bool same_as(const Wall & that) const;
//
//        virtual Vectorf<3> normal_reflect(
//                const Vectorf<3> &start,
//                const Vectorf<3> &hit_point,
//                float length) const;
//
//        virtual Vectorf<3> normal_reflect(const Vectorf<3> &incident) const;
//
//        virtual float cosine_angle(   // cosine angle with respect to surface normal
//                const Vectorf<3> &p
//        ) const;
};
//
class SimplePolygon : public Polygon
{
private:

    // Wall geometry properties
    Eigen::Matrix<float, 3, 1>  normal;
    Eigen::Matrix<float, 3, Eigen::Dynamic> corners;

    /* for 3D wall, provide local basis for plane of wall */
    Eigen::Matrix<float, 3, 1> origin;
    Eigen::Matrix<float, 3, 2> basis;
    Eigen::Matrix<float, 2, Eigen::Dynamic> flat_corners;

public:

    // Constructor
//    SimplePolygon(
//            const Eigen::Matrix<float, D, Eigen::Dynamic> &_corners,
//            const Eigen::ArrayXf &_absorption,
//            const Eigen::ArrayXf &_scatter,
//            const std::string &_name
//    );
//    Wall(
//            const Eigen::Matrix<float, D, Eigen::Dynamic> &_corners,
//            const Eigen::ArrayXf &_absorption,
//            const Eigen::ArrayXf &_scatter
//    ) : Wall(_corners, _absorption, _scatter, "") {}

    virtual float area() const;  // compute the area of the wall
    virtual Eigen::Matrix<float, 3, 1> get_origin() const;
    virtual int intersection(  // compute the intersection of line segment (p1 <-> p2) with wall
            const Vectorf<3> &p1,
            const Vectorf<3> &p2,
            Eigen::Ref<Vectorf<3>> intersection
    ) const;

//    virtual int reflect(
//            const Vectorf<D> &p,
//            Eigen::Ref<Vectorf<D>> p_reflected
//    ) const;
//    virtual int side(const Vectorf<D> &p) const;
//    virtual bool same_as(const Wall & that) const;
//
//    virtual Vectorf<D> normal_reflect(
//            const Vectorf<D> &start,
//            const Vectorf<D> &hit_point,
//            float length) const;
//
//    virtual Vectorf<D> normal_reflect(const Vectorf<D> &incident) const;
//
//    virtual float cosine_angle(   // cosine angle with respect to surface normal
//            const Vectorf<D> &p
//    ) const;

    SimplePolygon(
            Eigen::Matrix<float, 3, Eigen::Dynamic> &_corners,
            Eigen::Matrix<float, 3, 1> &_origin
    );

    // _corners.col(0) is guaranteed to be the correct type by the type constraint of _corners
    SimplePolygon(
            Eigen::Matrix<float, 3, Eigen::Dynamic> &_corners
    ) : SimplePolygon(_corners, (Eigen::Matrix<float, 3, 1>) _corners.col(0)) {}

};

class PolygonWithHole : public Polygon
{

private:
    SimplePolygon outer_polygon;
    std::vector<SimplePolygon> inner_polygons;

public:
    virtual float area() const;  // compute the area of the wall
    virtual Eigen::Matrix<float, 3, 1> get_origin() const;
    virtual int intersection(  // compute the intersection of line segment (p1 <-> p2) with wall
            const Vectorf<3> &p1,
            const Vectorf<3> &p2,
            Eigen::Ref<Vectorf<3>> intersection
    ) const;

    PolygonWithHole(
            Eigen::Matrix<float, 3, Eigen::Dynamic> &_corners,
            std::vector<Eigen::Matrix<float, 3, Eigen::Dynamic>> &_holes
    );
};

template<size_t D>
class Wall
{
  protected:
    void init();  // common part of initialization for walls of any dimension

  public:
    enum Isect {  // The different cases for intersections
      NONE = -1,  // - There is no intersection
      VALID = 0,  // - There is a valid intersection
      ENDPT = 1,  // - The intersection is on the endpoint of the segment
      BNDRY = 2   // - The intersection is on the boundary of the wall
    };

    static const int dim = D;

    // Wall properties container
    Eigen::ArrayXf absorption;  // the wall absorption coefficient for every freq. band
    Eigen::ArrayXf scatter;  // the wall scattering coefficient for every freq. band
    std::string name;
    Eigen::ArrayXf transmission;  // computed from absorption as sqrt(1 - a)
    Eigen::ArrayXf energy_reflection;  // computed from absorption as (1 - a)
    
    // Wall geometry properties
    Eigen::Matrix<float, D, 1>  normal;
    Eigen::Matrix<float, D, 1> origin;
    Eigen::Matrix<float, D, Eigen::Dynamic> corners;

    // Constructor
    Wall(
        const Eigen::Matrix<float,D,Eigen::Dynamic> &_corners,
        const Eigen::ArrayXf &_absorption,
        const Eigen::ArrayXf &_scatter,
        const std::string &_name
        );

//    // Copy constructor
//    Wall(const Wall<D> &w) :
//      absorption(w.absorption), scatter(w.scatter), name(w.name),
//      transmission(w.transmission), energy_reflection(w.energy_reflection),
//      normal(w.normal), corners(w.corners),
//      origin(w.origin), basis(w.basis), flat_corners(w.flat_corners)
//    {}

    // public methods

    // getters
    const Eigen::ArrayXf &get_transmission() const { return transmission; }
    const Eigen::ArrayXf &get_energy_reflection() const { return energy_reflection; }
    size_t get_n_bands() const { return transmission.size(); }

    // methods specific to the concrete types of Wall<D>
    virtual float area() const;  // compute the area of the wall
    virtual int intersection(  // compute the intersection of line segment (p1 <-> p2) with wall
        const Vectorf<D> &p1,
        const Vectorf<D> &p2,
        Eigen::Ref<Vectorf<D>> intersection
        ) const;

    // utility methods
    int intersects(
        const Vectorf<D> &p1,
        const Vectorf<D> &p2
        ) const;

    int reflect(
        const Vectorf<D> &p,
        Eigen::Ref<Vectorf<D>> p_reflected
        ) const;
    int side(const Vectorf<D> &p) const;
    bool same_as(const Wall & that) const;

    Vectorf<D> normal_reflect(
        const Vectorf<D> &start,
        const Vectorf<D> &hit_point,
        float length) const;

    Vectorf<D> normal_reflect(const Vectorf<D> &incident) const;

    float cosine_angle(   // cosine angle with respect to surface normal
        const Vectorf<D> &p
        ) const;
};

class Wall2D: public Wall<2>
{

    public:


        virtual float area() const;  // compute the area of the wall
        virtual int intersection(  // compute the intersection of line segment (p1 <-> p2) with wall
                const Vectorf<2> &p1,
                const Vectorf<2> &p2,
                Eigen::Ref<Vectorf<2>> intersection
        ) const;

        // Constructor
        Wall2D(
                const Eigen::Matrix<float, 2, Eigen::Dynamic> &_corners,
                const Eigen::ArrayXf &_absorption,
                const Eigen::ArrayXf &_scatter,
                const std::string &_name
        );
        Wall2D(
                const Eigen::Matrix<float, 2, Eigen::Dynamic> &_corners,
                const Eigen::ArrayXf &_absorption,
                const Eigen::ArrayXf &_scatter
        ) : Wall2D(_corners, _absorption, _scatter, "") {}

};

class Wall3D: public Wall<3>
{

    public:

        /* for 3D wall, provide local basis for plane of wall */
        Eigen::Matrix<float, 3, 2> basis;
        Eigen::Matrix<float, 2, Eigen::Dynamic> flat_corners;
//        std::unique_ptr<Polygon> wall_geometry;
        Polygon* wall_geometry;

        virtual float area() const;  // compute the area of the wall
        virtual int intersection(  // compute the intersection of line segment (p1 <-> p2) with wall
                const Vectorf<3> &p1,
                const Vectorf<3> &p2,
                Eigen::Ref<Vectorf<3>> intersection
        ) const;

        // Constructor
        Wall3D(
                const Eigen::Matrix<float, 3, Eigen::Dynamic> &_corners,
                const std::vector<Eigen::Matrix<float, 3, Eigen::Dynamic>> &_holes,
                const Eigen::ArrayXf &_absorption,
                const Eigen::ArrayXf &_scatter,
                const std::string &_name
        );
        Wall3D(
                const Eigen::Matrix<float, 3, Eigen::Dynamic> &_corners,
                const std::vector<Eigen::Matrix<float, 3, Eigen::Dynamic>> &_holes,
                const Eigen::ArrayXf &_absorption,
                const Eigen::ArrayXf &_scatter
        ) : Wall3D(_corners, _holes, _absorption, _scatter, "") {}
        Wall3D(
                const Eigen::Matrix<float, 3, Eigen::Dynamic> &_corners,
                const Eigen::ArrayXf &_absorption,
                const Eigen::ArrayXf &_scatter,
                const std::string &_name
        ) : Wall3D(
                _corners,
                std::vector<Eigen::Matrix<float, 3, Eigen::Dynamic>>(),
                _absorption,
                _scatter,
                _name){}
        Wall3D(
                const Eigen::Matrix<float, 3, Eigen::Dynamic> &_corners,
                const Eigen::ArrayXf &_absorption,
                const Eigen::ArrayXf &_scatter
        ) : Wall3D(
                _corners,
                std::vector<Eigen::Matrix<float, 3, Eigen::Dynamic>>(),
                _absorption,
                _scatter,
                ""){}

        ~Wall3D()
        {
            delete wall_geometry;
        }
};


#include "wall.cpp"

#endif // __CWALL_H__
