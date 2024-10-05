#include <gtest/gtest.h>

#include "mytypes.h"
#include "vector.c"


#define VEC_SIZE (100)


namespace
{
    class VectorTest : public ::testing::Test
    {
        protected:
            real *a;
            real *b;

            VectorTest( )
            {
                if ( (a = (real *) malloc( VEC_SIZE * sizeof(real))) == NULL ||
                            (b = (real *) malloc( VEC_SIZE * sizeof(real))) == NULL )
                {
                    throw new std::bad_alloc( );
                }
            }

            virtual ~VectorTest( )
            {
                if ( a != NULL )
                {
                    free( a );
                }
                if ( b != NULL )
                {
                    free( b );
                }
            }

            virtual void SetUp( )
            {
                for ( int i = 0; i < VEC_SIZE; ++i )
                {
                    a[i] = i + 1.0;
                    b[i] = 1.0;
                }
            }

            virtual void TearDown( )
            {

            }
    };


    TEST_F(VectorTest, Dot)
    {
        ASSERT_EQ( Dot(a, b, VEC_SIZE), (VEC_SIZE * (VEC_SIZE + 1.0)) / 2.0 );
        ASSERT_EQ( Dot(b, b, VEC_SIZE), (real) VEC_SIZE );
    }


    TEST_F(VectorTest, Norm)
    {
        ASSERT_EQ( Norm(a, VEC_SIZE), sqrt( VEC_SIZE * (VEC_SIZE + 1.0) * (2.0 * VEC_SIZE + 1.0) / 6.0 ) );
        ASSERT_EQ( Norm(b, VEC_SIZE), sqrt( (real) VEC_SIZE ) );
    }
}


int main( int argc, char **argv )
{
    ::testing::InitGoogleTest( &argc, argv );
    return RUN_ALL_TESTS( );
}
