#include <gtest/gtest.h>

#include "mytypes.h"
#include "allocate.c"
#include "lin_alg.c"
#include "list.c"
#include "tool_box.c"
#include "vector.c"

#define VEC_SIZE (100)


namespace
{
    class LinAlgTest : public ::testing::Test
    {
        protected:
            real *a, *res;
            sparse_matrix *sm;

            LinAlgTest( )
            {
                if ( !Allocate_Matrix( &sm, VEC_SIZE, VEC_SIZE ) ||
                        (a = (real *) malloc( VEC_SIZE * sizeof(real))) == NULL ||
                        (res = (real *) malloc( VEC_SIZE * sizeof(real))) == NULL )
                {
                    throw new std::bad_alloc( );
                }
            }

            virtual ~LinAlgTest( )
            {
                if ( a != NULL )
                {
                    free( a );
                }
                if ( res != NULL )
                {
                    free( res );
                }
                if ( sm != NULL )
                {
                    Deallocate_Matrix( sm );
                }

            }

            virtual void SetUp( )
            {
                for ( int i = 0; i < VEC_SIZE; ++i )
                {
                    a[i] = 1.0;
                }

                // set up sparse matrix which is an identity matrix in our case
                for ( int i = 0; i < sm->n + 1; i++ )
                {
                    sm->start[i] = i;
                }

                for ( int i = 0; i < sm->start[sm->n]; i++ )
                {
                    sm->j[i] = i;
                    sm->val[i] = 1.0;
                }
            }


            virtual void TearDown( )
            {

            }
    };


    TEST_F(LinAlgTest, Sparse_MatVec)
    {
        Sparse_MatVec( sm, a, res );

        for ( int i = 0; i < VEC_SIZE; ++i )
        {
            ASSERT_EQ( res[i], a[i] );
        }
    }
}


int main( int argc, char **argv )
{
    ::testing::InitGoogleTest( &argc, argv );
    return RUN_ALL_TESTS( );
}
