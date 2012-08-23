//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Copyright 2012 Matthew R. Renaud
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//m=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This file is part of Protein-protein-interaction.
//
// Protein-protein-interaction is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation, either
// version 3 of the License, or(at your option) any later version.
//
// Protein-protein-interaction is distributed in the hope that it will
// be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Protein-protein-interaction. If not, see
// <http://www.gnu.org/licenses/>.
// =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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
