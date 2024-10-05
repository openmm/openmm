#include <gtest/gtest.h>

#include "spuremd.h"


namespace
{
    class SPuReMDTest : public ::testing::Test
    {
        protected:
            void *handle;

            SPuReMDTest ( )
            {
            }

            virtual ~SPuReMDTest ( )
            {
            }

            virtual void SetUp( )
            {
            }

            virtual void TearDown( )
            {
                if ( handle != NULL )
                {
                    cleanup( handle );
                }
            }
    };


    TEST_F(SPuReMDTest, water_6540)
    {
        handle = setup( "../data/benchmarks/water/water_6540.pdb", 
                "../data/benchmarks/water/ffield.water",
                "../environ/param.gpu.water" );

        ASSERT_EQ( simulate( handle ), SPUREMD_SUCCESS );

        //TODO: check energy after evolving system, e.g., 100 steps
    }


    TEST_F(SPuReMDTest, silica_6000)
    {
        handle = setup( "../data/benchmarks/silica/silica_6000.pdb", 
                "../data/benchmarks/silica/ffield-bio",
                "../environ/param.gpu.water" );

        ASSERT_EQ( simulate( handle ), SPUREMD_SUCCESS );

        //TODO: check energy after evolving system, e.g., 100 steps
    }
}


int main( int argc, char **argv )
{
    ::testing::InitGoogleTest( &argc, argv );
    return RUN_ALL_TESTS( );
}
