#include "basis.hpp"
#include "lattice.hpp"
#include "solver.hpp"
#include "spinOPft.hpp"
#include <gtest/gtest.h>

/* ====================================================== */
/* NOTE: The global variable approach for argc and argv is 
   VERY important to pass the command line arguments to the 
   test (like -malloc_dump). */
/* ====================================================== */

static char help[] = "This is a test to test the Sqp operator class\n\n";
int _argc;
char ** _argv;

Environment env{_argc,_argv,6,help};

class HamiltonianTestEnv : public ::testing::Environment {
protected:

  virtual void TearDown() {
    ::testing::internal::CaptureStdout();
    EXPECT_EQ("", ::testing::internal::GetCapturedStdout());
  }
};


TEST(SqpOperator, chain1D) {

  Basis b0{env,0};
  Basis b1{env,1};
  // for (int i = 0; i < b.local_size; ++i)
  //   PetscSynchronizedPrintf(PETSC_COMM_WORLD, "basis[%d]: %d\n", i, b.int_basis[i]);
  // PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

  PetscErrorCode ierr = 0;

  Sq_data sqp_data1{&b0,&b1};
  Sqp sqp1{env, sqp_data1};
  
  if (b0.mpi_rank == 0) {
    EXPECT_EQ(7, b0.int_basis[0]);
    sqp1.OpOnBasisElems(env,b0.int_basis[0],0);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[0],1);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[0],2);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[0],3);
    EXPECT_NEAR(0, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[0],4);
    EXPECT_NEAR(1, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[0],5);
    EXPECT_NEAR(5, sqp1.value, 1e-13);

    EXPECT_EQ(11, b0.int_basis[1]);
    sqp1.OpOnBasisElems(env,b0.int_basis[1],0);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[1],1);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[1],2);
    EXPECT_NEAR(0, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[1],3);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[1],4);
    EXPECT_NEAR(2, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[1],5);
    EXPECT_NEAR(6, sqp1.value, 1e-13);
  }

  if (b0.mpi_rank == 1) {
    EXPECT_EQ(35, b0.int_basis[3]);
    sqp1.OpOnBasisElems(env,b0.int_basis[3],0);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[3],1);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[3],2);
    EXPECT_NEAR(5, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[3],3);
    EXPECT_NEAR(6, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[3],4);
    EXPECT_NEAR(9, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[3],5);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);

    EXPECT_EQ(41, b0.int_basis[6]);
    sqp1.OpOnBasisElems(env,b0.int_basis[6],0);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[6],1);
    EXPECT_NEAR(6, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[6],2);
    EXPECT_NEAR(7, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[6],3);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[6],4);
    EXPECT_NEAR(12, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[6],5);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);
  }

  if (b0.mpi_rank == 2) {
    EXPECT_EQ(49, b0.int_basis[2]);
    sqp1.OpOnBasisElems(env,b0.int_basis[2],0);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[2],1);
    EXPECT_NEAR(9, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[2],2);
    EXPECT_NEAR(10, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[2],3);
    EXPECT_NEAR(12, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[2],4);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[2],5);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);

    EXPECT_EQ(52, b0.int_basis[4]);
    sqp1.OpOnBasisElems(env,b0.int_basis[4],0);
    EXPECT_NEAR(10, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[4],1);
    EXPECT_NEAR(11, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[4],2);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[4],3);
    EXPECT_NEAR(14, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[4],4);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);
    sqp1.OpOnBasisElems(env,b0.int_basis[4],5);
    EXPECT_NEAR(-1, sqp1.value, 1e-13);
    
  }

  PetscInt * index1;
  PetscScalar * values1;
  ierr = PetscCalloc1(5, &index1); ASSERT_EQ(0,ierr);
  ierr = PetscCalloc1(5, &values1); ASSERT_EQ(0,ierr);
  index1[0] = 4;
  index1[1] = 6;
  index1[2] = 11;
  index1[3] = 15;
  index1[4] = 19;

  for (PetscInt i = 0; i < 5; ++i)
    values1[i] = 2.0;

  Vec state1, rhs1;
  const PetscScalar * rhs_array1;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, b0.size, &state1); ASSERT_EQ(0,ierr);
  ierr = VecSetValues(state1, 5, index1, values1, INSERT_VALUES); ASSERT_EQ(0,ierr);
  ierr = VecAssemblyBegin(state1); ASSERT_EQ(0,ierr);
  ierr = VecAssemblyEnd(state1); ASSERT_EQ(0,ierr);
  ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, b1.size, &rhs1); ASSERT_EQ(0,ierr);

  ierr = sqp1.OpOnStateVector(env, state1, rhs1, 0); ASSERT_EQ(0,ierr);
  ierr = VecGetArrayRead(rhs1, &rhs_array1); ASSERT_EQ(0,ierr);
  
  if (b1.mpi_rank == 0) {
    EXPECT_NEAR(0.0, PetscRealPart(rhs_array1[0]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array1[0]), 1e-13);

    EXPECT_NEAR(2.0, PetscRealPart(rhs_array1[1]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array1[1]), 1e-13);

    EXPECT_NEAR(0.0, PetscRealPart(rhs_array1[2]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array1[2]), 1e-13);
    
    EXPECT_NEAR(0.0, PetscRealPart(rhs_array1[3]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array1[3]), 1e-13);
    
    EXPECT_NEAR(0.0, PetscRealPart(rhs_array1[4]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array1[4]), 1e-13);
  }

  if (b1.mpi_rank == 1) {
    EXPECT_NEAR(0.0, PetscRealPart(rhs_array1[0]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array1[0]), 1e-13);

    EXPECT_NEAR(0.0, PetscRealPart(rhs_array1[1]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array1[1]), 1e-13);

    EXPECT_NEAR(2.0, PetscRealPart(rhs_array1[2]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array1[2]), 1e-13);
    
    EXPECT_NEAR(0.0, PetscRealPart(rhs_array1[3]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array1[3]), 1e-13);
    
    EXPECT_NEAR(0.0, PetscRealPart(rhs_array1[4]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array1[4]), 1e-13);
  }

  if (b1.mpi_rank == 2) {
    EXPECT_NEAR(0.0, PetscRealPart(rhs_array1[0]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array1[0]), 1e-13);

    EXPECT_NEAR(0.0, PetscRealPart(rhs_array1[1]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array1[1]), 1e-13);

    EXPECT_NEAR(2.0, PetscRealPart(rhs_array1[2]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array1[2]), 1e-13);
    
    EXPECT_NEAR(0.0, PetscRealPart(rhs_array1[3]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array1[3]), 1e-13);
    
    EXPECT_NEAR(0.0, PetscRealPart(rhs_array1[4]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array1[4]), 1e-13);
  }
 
  ierr = VecRestoreArrayRead(rhs1, &rhs_array1); ASSERT_EQ(0,ierr);
  ierr = PetscFree(index1); ASSERT_EQ(0,ierr);
  ierr = PetscFree(values1); ASSERT_EQ(0,ierr);
  ierr = PetscFree(rhs_array1); ASSERT_EQ(0,ierr);
  ierr = VecDestroy(&state1); ASSERT_EQ(0,ierr);
  ierr = VecDestroy(&rhs1); ASSERT_EQ(0,ierr);

  // ----------------------------------------------------------------------------------------- //
  // ----------------------------------------------------------------------------------------- //

  Sq_data sqp_data2{&b0,&b1};
  Sqp sqp2{env, sqp_data2};
  
  PetscInt * index2;
  PetscScalar * values2;
  ierr = PetscCalloc1(5, &index2); ASSERT_EQ(0,ierr);
  ierr = PetscCalloc1(5, &values2); ASSERT_EQ(0,ierr);
  index2[0] = 0;
  index2[1] = 2;
  index2[2] = 9;
  index2[3] = 13;
  index2[4] = 17;

  for (PetscInt i = 0; i < 5; ++i)
    values2[i] = 0.25;

  Vec state2, rhs2;
  const PetscScalar * rhs_array2;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, b0.size, &state2); ASSERT_EQ(0,ierr);
  ierr = VecSetValues(state2, 5, index2, values2, INSERT_VALUES); ASSERT_EQ(0,ierr);
  ierr = VecAssemblyBegin(state2); ASSERT_EQ(0,ierr);
  ierr = VecAssemblyEnd(state2); ASSERT_EQ(0,ierr);
  ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, b1.size, &rhs2); ASSERT_EQ(0,ierr);
  
  ierr = sqp1.OpOnStateVector(env, state2, rhs2, 3); ASSERT_EQ(0,ierr);  
  ierr = VecGetArrayRead(rhs2, &rhs_array2); ASSERT_EQ(0,ierr);
  
  if (b1.mpi_rank == 0) {
    EXPECT_NEAR(0.25, PetscRealPart(rhs_array2[0]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array2[0]), 1e-13);

    EXPECT_NEAR(0.0, PetscRealPart(rhs_array2[1]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array2[1]), 1e-13);

    EXPECT_NEAR(0.0, PetscRealPart(rhs_array2[2]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array2[2]), 1e-13);
    
    EXPECT_NEAR(0.0, PetscRealPart(rhs_array2[3]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array2[3]), 1e-13);
    
    EXPECT_NEAR(0.0, PetscRealPart(rhs_array2[4]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array2[4]), 1e-13);
  }

  if (b1.mpi_rank == 1) {
    EXPECT_NEAR(0.0, PetscRealPart(rhs_array2[0]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array2[0]), 1e-13);

    EXPECT_NEAR(0.0, PetscRealPart(rhs_array2[1]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array2[1]), 1e-13);

    EXPECT_NEAR(0.0, PetscRealPart(rhs_array2[2]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array2[2]), 1e-13);
    
    EXPECT_NEAR(0.0, PetscRealPart(rhs_array2[3]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array2[3]), 1e-13);
    
    EXPECT_NEAR(0.0, PetscRealPart(rhs_array2[4]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array2[4]), 1e-13);
  }

  if (b1.mpi_rank == 2) {
    EXPECT_NEAR(0.0, PetscRealPart(rhs_array2[0]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array2[0]), 1e-13);

    EXPECT_NEAR(0.0, PetscRealPart(rhs_array2[1]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array2[1]), 1e-13);

    EXPECT_NEAR(0.0, PetscRealPart(rhs_array2[2]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array2[2]), 1e-13);
    
    EXPECT_NEAR(0.25, PetscRealPart(rhs_array2[3]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array2[3]), 1e-13);
    
    EXPECT_NEAR(0.0, PetscRealPart(rhs_array2[4]), 1e-13);
    EXPECT_NEAR(0.0, PetscImaginaryPart(rhs_array2[4]), 1e-13);
  }
 
  ierr = VecRestoreArrayRead(rhs2, &rhs_array2); ASSERT_EQ(0,ierr);
  ierr = PetscFree(index2); ASSERT_EQ(0,ierr);
  ierr = PetscFree(values2); ASSERT_EQ(0,ierr);
  ierr = PetscFree(rhs_array2); ASSERT_EQ(0,ierr);
  ierr = VecDestroy(&state2); ASSERT_EQ(0,ierr);
  ierr = VecDestroy(&rhs2); ASSERT_EQ(0,ierr);
}


int main (int argc, char ** argv){
  _argc = argc;
  _argv = argv;
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(new HamiltonianTestEnv);
  return RUN_ALL_TESTS();
}

