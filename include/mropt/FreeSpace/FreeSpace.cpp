#include "FreeSpace.hpp"  

using namespace mropt::freespace;

FreeSpace::~FreeSpace()
{
	if(poly_centers) delete poly_centers;
}

int FreeSpace::Polygon::counter = 0;
std::vector<FreeSpace::Polygon> FreeSpace::polygons = std::vector<FreeSpace::Polygon>();
std::unordered_map<
    std::pair<int, int>,
        mropt::StateSpace::State::state,
        boost::hash<std::pair<int, int>>>*
FreeSpace::poly_centers = nullptr;


void FreeSpace::add_polygon(MX &A, MX &b)
{
  polygons.push_back(Polygon{A, b});
}

void FreeSpace::add(std::list<Polygon> polys)
{
  std::copy(polys.begin(), polys.end(), std::back_inserter(polygons));
}


std::vector<MX> FreeSpace::get_constraints(std::vector<PolygonAssignment> pas)
{
  double safe_radius = robot_shape->get_safety_radius();
  MX zero{0.0};
  // 用于累积所有约束条件的总和
  MX g_sum{0.0};
  std::vector<MX> constraints{};
  // 获得状态变量
  const auto &xy = ss.xy();
  for (const auto &pa : pas)
  {
    // 获得多边形的信息
    auto &polygon = polygons[pa.pid];
    for(int k = pa.k0; k < pa.kf; ++k){
      // mtimes是符号计算, g_p存储符号计算结果, 这里的变量就是xy
      auto g_p = mtimes(polygon.A, xy(all, k)) - polygon.b + safe_radius + threshold;
      // 将计算结果加入向量
      constraints.push_back( g_p );
      for (int g_id = 0; g_id < g_p.size1(); ++g_id)
      {
        g_sum = g_sum + MX::mmax(MX::vertcat({g_p(g_id), zero}));
      }
    }
  }
  // 构成所有约束条件的累加和, 传入参数能够获得该累加和的值
  J_real_ = Function("J_real", {ss.X()}, {g_sum});
  return constraints;
}

void FreeSpace::clear_polygons() {
  FreeSpace::Polygon::counter = 0;
  FreeSpace::polygons = std::vector<FreeSpace::Polygon>();
  delete FreeSpace::poly_centers;
  FreeSpace::poly_centers = nullptr;
}