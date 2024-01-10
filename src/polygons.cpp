#include "../include/pendulum_feasibility_solver/polygons.h"
// clang-format off
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>
#include <boost/geometry/multi/geometries/register/multi_polygon.hpp>
// clang-format on

BOOST_GEOMETRY_REGISTER_BOOST_TUPLE_CS(cs::cartesian)
namespace bg = boost::geometry;

void Polygon::convex_hull()
{
  Polygone_Corners.clear();

  typedef boost::tuple<double, double> point;
  typedef bg::model::multi_point<point> points;
  typedef bg::model::polygon<point> polygon;

  points p;
  for(size_t r = 0; r < _Rectangles.size(); r++)
  {
    std::vector<Eigen::Vector3d> rect_corners = _Rectangles[r].Get_corners();
    for(size_t i = 0; i < rect_corners.size(); i++) { bg::append(p, point(rect_corners[i].x(), rect_corners[i].y())); }
  }
  polygon hull;
  bg::convex_hull(p, hull);

  for(auto it = boost::begin(boost::geometry::exterior_ring(hull));
      it != boost::end(boost::geometry::exterior_ring(hull)); ++it)
  {
    Polygone_Corners.push_back(Eigen::Vector3d{bg::get<0>(*it), bg::get<1>(*it), 0.});
  }
}

void Polygon::jarvis_march()
{
  // std::sort(_corners.begin(),_corners.end(),vec3d_x_comp());

  size_t index = std::min_element(_corners.begin(), _corners.end(), vec3d_x_comp()) - _corners.begin();

  Polygone_Corners.push_back(_corners[index]);

  size_t l = index;
  size_t q = 0;
  while(true)
  {
    q = (l + 1) % _corners.size();
    for(size_t i = 0; i < _corners.size(); i++)
    {
      if(i == l) { continue; }
      const Eigen::Vector3d & v_q{_corners[q].x(), _corners[q].y(), 0};
      const Eigen::Vector3d & v_l{_corners[l].x(), _corners[l].y(), 0};
      const Eigen::Vector3d & v_i{_corners[i].x(), _corners[i].y(), 0};
      double d = ((v_q - v_l).cross(v_i - v_l)).z();
      if(d > 0 || (d == 0 && (v_i - v_l).norm() > (v_q - v_l).norm())) { q = i; }
    }
    l = q;
    if(l == index) { break; }
    Polygone_Corners.push_back(_corners[q]);
  }
}
