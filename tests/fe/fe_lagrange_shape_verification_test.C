// libmesh includes
#include "libmesh/libmesh.h"
#include "libmesh/fe.h"

// unit test includes
#include "test_comm.h"

#include "libmesh_cppunit.h"

using namespace libMesh;

/**
 * This class is for unit testing lagrange shape function values and derivatives
 */
class LagrangeShapeTest : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE( LagrangeShapeTest );
  CPPUNIT_TEST( test1DLagrange );
  CPPUNIT_TEST_SUITE_END();

  void test1DLagrange ()
  {
    LOG_UNIT_TEST;

    // TOLERANCE*TOLERANCE works with double but not float128
    Real my_tol = TOLERANCE*std::sqrt(TOLERANCE);

    for (auto x : {0., 1./3., -.5, -1./7.})
      for(auto o : make_range(1, 4))
        for(auto i : make_range(o + 1))
        {
          const unsigned ii = (i > 0) ? i % o + 1: 0;
          LIBMESH_ASSERT_FP_EQUAL(fe_lagrange_1D_shape(Order(o), ii, x),
                                  fe_lagrange_1D_arbitrary_shape(Order(o), i, x), my_tol);
          LIBMESH_ASSERT_FP_EQUAL(fe_lagrange_1D_shape_deriv(Order(o), ii, 0, x),
                                  fe_lagrange_1D_arbitrary_shape_deriv(Order(o), i, 0, x), my_tol);
          LIBMESH_ASSERT_FP_EQUAL(fe_lagrange_1D_shape_second_deriv(Order(o), ii, 0, x),
                                  fe_lagrange_1D_arbitrary_shape_second_deriv(Order(o), i, 0, x), my_tol);
        }
  }

  void setUp() {}

  void tearDown() {}

};

CPPUNIT_TEST_SUITE_REGISTRATION( LagrangeShapeTest );
