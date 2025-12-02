#include <mcpp_heisenberg/geometry.hpp>

#include "tests.hpp"

class GeometryTestsSuite : public MCHTestsSuite {
 protected:
  GeometryTestsSuite() = default;
};

TEST_F(GeometryTestsSuite, TestCreate) {
  auto geometry = mch::Geometry();
}
