#pragma once
#include <Eigen/StdVector>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <SpaceVecAlg/SpaceVecAlg>

inline Eigen::Vector3d rpyFromMat(const Eigen::Matrix3d & E)
{
  double roll = atan2(E(1, 2), E(2, 2));
  double pitch = -asin(E(0, 2));
  double yaw = atan2(E(0, 1), E(0, 0));
  return Eigen::Vector3d(roll, pitch, yaw);
}

struct Rectangle
{

public:
  Rectangle(double ori, const Eigen::Vector2d size, const Eigen::Vector3d offset = Eigen::Vector3d::Zero())
  {
    _center = offset;
    _angle = ori;
    _size.segment(0, 2) = size;
    compute_rect();
  }

  Rectangle(const sva::PTransformd pose,
            const Eigen::Vector2d size,
            const Eigen::Vector3d offset = Eigen::Vector3d::Zero())
  {
    _center = pose.translation() + offset;
    _center.z() = 0;
    _angle = rpyFromMat(pose.rotation()).z();
    _size.segment(0, 2) = size;
    compute_rect();
  }

  Rectangle(const Eigen::Vector3d center,
            const Eigen::Vector2d size,
            const Eigen::Vector3d offset = Eigen::Vector3d::Zero())
  {
    _center = center + offset;
    _angle = _center.z();
    _size.segment(0, 2) = size;
    _center.z() = 0;
    compute_rect();
  }
  void compute_rect()
  {

    R << cos(_angle), -sin(_angle), 0, sin(_angle), cos(_angle), 0, 0, 0, 1;

    upper_left_corner = _center + R * Eigen::Vector3d{-_size.x() / 2, _size.y() / 2, 0};
    upper_right_corner = _center + R * Eigen::Vector3d{_size.x() / 2, _size.y() / 2, 0};
    lower_left_corner = _center + R * Eigen::Vector3d{-_size.x() / 2, -_size.y() / 2, 0};
    lower_right_corner = _center + R * Eigen::Vector3d{_size.x() / 2, -_size.y() / 2, 0};
    corners = {upper_left_corner, upper_right_corner, lower_right_corner, lower_left_corner};
  }
  void add_offset(const Eigen::Vector3d offset)
  {
    _center += Eigen::Vector3d{offset.x(), offset.y(), 0};
    compute_rect();
  }
  ~Rectangle()
  {
    corners.clear();
    _center.setZero();
    _size.setZero();
    _angle = 0;
    R = Eigen::Matrix3d::Identity();
    upper_left_corner.setZero();
    upper_right_corner.setZero();
    lower_left_corner.setZero();
    lower_right_corner.setZero();
  }

  std::vector<Eigen::Vector3d> & Get_corners()
  {
    return corners;
    // return {upper_left_corner,lower_left_corner,lower_right_corner,upper_right_corner};
  }
  const Eigen::Vector3d & Up_Left_corner() const noexcept
  {
    return upper_left_corner;
  }
  const Eigen::Vector3d & Up_Right_corner() const noexcept
  {
    return upper_right_corner;
  }
  const Eigen::Vector3d & Dwn_Right_corner() const noexcept
  {
    return lower_right_corner;
  }
  const Eigen::Vector3d & Dwn_Left_corner() const noexcept
  {
    return lower_left_corner;
  }

  double get_yaw()
  {
    return _angle;
  }
  Eigen::Vector3d get_center()
  {
    return _center;
  }

private:

  Eigen::Vector3d _center = Eigen::Vector3d::Zero();
  Eigen::Vector3d _size = Eigen::Vector3d::Zero();
  double _angle = 0;
  Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
  Eigen::Vector3d upper_left_corner = Eigen::Vector3d::Zero();
  Eigen::Vector3d upper_right_corner = Eigen::Vector3d::Zero();
  Eigen::Vector3d lower_left_corner = Eigen::Vector3d::Zero();
  Eigen::Vector3d lower_right_corner = Eigen::Vector3d::Zero();
  std::vector<Eigen::Vector3d> corners;
};

struct vec3d_x_comp
{
  inline bool operator()(const Eigen::Vector3d & struct1, const Eigen::Vector3d & struct2)
  {
    return (struct1.x() < struct2.x());
  }
};

struct Polygon
{

public:
  Polygon() = default;
  Polygon(const Rectangle Rect1, const Rectangle Rect2)
  {
    _Rectangles = {Rect1, Rect2};
    Compute_polygone();
  }
  Polygon(const Rectangle Rect1)
  {
    _Rectangles = {Rect1};
    Compute_polygone();
  }
  Polygon(const Eigen::MatrixX2d normals, Eigen::VectorXd offsets)
  {
    Polygone_Normals = normals;
    Offset = offsets;
    cstr_to_polygone();
  }

  ~Polygon()
  {
    _Rectangles.clear();
    _corners.clear();
    Polygone_Corners.clear();
  }

  void jarvis_march();
  void convex_hull();

  std::vector<Eigen::Vector3d> & Get_Polygone_Corners()
  {
    return Polygone_Corners;
  }

  const Eigen::MatrixX2d & normals()
  {
    return Polygone_Normals;
  }

  const Eigen::VectorXd & offsets()
  {
    return Offset;
  }

  Eigen::Vector3d get_center()
  {
    Eigen::Vector3d center = Eigen::Vector3d::Zero();
    double n = static_cast<double>(_Rectangles.size());
    for (auto & r : _Rectangles)
    {
      center += (1/n) * r.get_center();
    }
    return center;
  }

  Rectangle & get_Rectangle(int indx)
  {
    if(indx < static_cast<int>(_Rectangles.size()) - 1)
    {
      return _Rectangles[indx];
    }
    else
    {
      return _Rectangles[0];
    }
  }

private:
  void Compute_polygone()
  {
    if(_Rectangles.size() > 1)
    {
      // for (int r = 0 ; r < _Rectangles.size() ; r++)
      // {
      //     std::vector<Eigen::Vector3d> corners = _Rectangles[r].Get_corners();
      //     _corners.insert(_corners.end(),corners.begin(),corners.end());
      // }
      // jarvis_march();
      convex_hull();
    }
    else
    {
      Polygone_Corners = _Rectangles[0].Get_corners();
    }
    Polygone_Normals.resize(Polygone_Corners.size(), 2);
    Polygone_Edges_Center.resize(Polygone_Corners.size(), 2);
    Polygone_Vertices.resize(Polygone_Corners.size(), 2);
    Offset.resize(Polygone_Corners.size());
    Eigen::Matrix2d R_Vertices_0;

    for(size_t c = 0; c < Polygone_Corners.size(); c++)
    {
      const Eigen::Vector3d & point_1 = Polygone_Corners[c];
      const Eigen::Vector3d & point_2 = Polygone_Corners[(c + 1) % Polygone_Corners.size()];
      const Eigen::Vector3d vertice = (point_2 - point_1).normalized();
      const Eigen::Vector3d normal = vertical_vec.cross(vertice);
      Polygone_Normals(c, 0) = normal.x();
      Polygone_Normals(c, 1) = normal.y();
      Polygone_Vertices(c, 0) = vertice.x();
      Polygone_Vertices(c, 1) = vertice.y();
      Polygone_Edges_Center(c, 0) = (((point_2 + point_1) / 2)).x();
      Polygone_Edges_Center(c, 1) = (((point_2 + point_1) / 2)).y();

      R_Vertices_0 << Polygone_Normals(c, 0), Polygone_Vertices(c, 0), Polygone_Normals(c, 1),
          Polygone_Vertices(c, 1);

      Offset(c) = (R_Vertices_0.transpose() * Polygone_Edges_Center.block(c, 0, 1, 2).transpose())(0);
    }
  }
  
  void cstr_to_polygone()
  {

    std::vector<Eigen::Index> vertices_indx;
    for(Eigen::Index i = 0 ; i < Polygone_Normals.rows() ; i++)
    {
      Eigen::Index end_indx = (i + 1) % static_cast<Eigen::Index>(Polygone_Normals.rows());
      Eigen::RowVector2d ni = Polygone_Normals.block(i,0,1,2).normalized();
      Eigen::RowVector2d nip1 = Polygone_Normals.block(end_indx,0,1,2).normalized();
      // mc_rtc::log::info("normal {}\n{}\nnext_normal{}\ndot prod {}",i,ni,nip1,ni.transpose() * nip1);
      if( std::abs( ni * nip1.transpose() - 1) > 1e-4  )
      {

        vertices_indx.push_back(i);
        // mc_rtc::log::info("selected");
      }
    }
    // mc_rtc::log::info("corner {} selected {}",Polygone_Normals.rows(),normals.size());
    for(size_t i = 0 ; i < vertices_indx.size() ; i++)
    {
      Eigen::Index start_indx = vertices_indx[i];
      Eigen::Index end_indx = vertices_indx[(i + 1) % static_cast<Eigen::Index>(vertices_indx.size())];
      Eigen::Matrix2d R = Eigen::Matrix2d::Zero();
      Eigen::Vector2d o = Eigen::Vector2d::Zero(); 
      R.block(0,0,1,2) = Polygone_Normals.block(start_indx,0,1,2);
      R.block(1,0,1,2) = Polygone_Normals.block(end_indx,0,1,2);
      o.x() = Offset[start_indx];
      o.y() = Offset[end_indx];
      if(R.determinant() != 0)
      {
        Eigen::Vector2d p = R.inverse() * o;
        // mc_rtc::log::info("start {} ; end {} point {}",i,end_indx,p);
        // mc_rtc::log::info("R {}\no {}",R,o);
  
        Polygone_Corners.push_back(Eigen::Vector3d{p.x(),p.y(),0});
      }
      
      
    
    }
  }

  Eigen::Vector3d vertical_vec = Eigen::Vector3d{0., 0., 1.};
  std::vector<Rectangle> _Rectangles;
  std::vector<Eigen::Vector3d> _corners;
  std::vector<Eigen::Vector3d> Polygone_Corners;
  Eigen::MatrixX2d Polygone_Vertices;
  Eigen::MatrixX2d Polygone_Edges_Center;
  Eigen::MatrixX2d Polygone_Normals;
  Eigen::VectorXd Offset;
};