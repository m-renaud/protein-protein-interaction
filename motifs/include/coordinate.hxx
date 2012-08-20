#ifndef MRR_COORDINATE_HXX_
#define MRR_COORDINATE_HXX_

#include <cmath>

namespace mrr {

template <typename PointType>
struct coordinate;

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
template <typename PointType = double>
struct coordinate
{
  using value_type = PointType;
  
  coordinate()
    : coordinate(0.0, 0.0, 0.0)
  {
  }
  
  coordinate(value_type x_, value_type y_, value_type z_)
    : x(x_), y(y_), z(z_)
  {
  }

  value_type x, y, z;
};


template <typename PointType>
PointType distance(
  mrr::coordinate<PointType> const& a,
  mrr::coordinate<PointType> const& b
)
{
  return std::sqrt(
    std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2) + std::pow(a.z - b.z, 2)
  );
}



} // namespace mrr


//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
template <typename T>
std::ostream& operator <<(std::ostream& os, mrr::coordinate<T> const& c)
{
  os << "(" << c.x << ", " << c.y << ", " << c.z << ")";
  return os;
}



#endif // #ifndef MRR_COORDINATE_HXX_
