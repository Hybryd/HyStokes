#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>


typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

using Eigen::MatrixXf;
using Eigen::Array44f;
using Eigen::Array22f;
using Eigen::VectorXi;
using Eigen::VectorXf;
using Eigen::VectorXd;

using namespace std;

/*
  Stokes equation for incompressible fluid:
  -----------------------------------------
  
  On the square:
  -\nu \Delta u + \nabla p = f
                    div u = 0
  On the boundary:
                        u = g

  We use here the finite difference method with euclidean grid to solve the problem.
  The grid:
  
  (N+1,0) - (N+1,1) - ... - (N+1,N+1)
      .         .               .
      .         .               .
      .         .               .
    (1,0) -   (1,1) - ... -   (1,N+1)
    (0,0) -   (0,1) - ... -   (0,N+1)
  
  The points are indexed the following way from the bottom left to the upper right, from 0 to (N+1)^2-1
  
  (N+1)(N+2) (N+1)(N+2)+1 ... (N+1)^2-1
      .           .               .
      .           .               .
      .           .               .
     N+2         N+3      ...  N+2+N+1
      0           1       ...    N+1
  
  Global indices  : from 0 to (N+1)^2-1
  Boundary indices: from 0 to 4*(N-1)
  Inner indices   : from 0 to N^2
  
  Need the Eigen library: http://eigen.tuxfamily.org/
  
  Compile with (-std=gnu++11 because of the Lambda calculus inside of prune)
    g++ *.cpp -o main -std=gnu++11

  Usage :
    ./main N > out
    where N is the number of points for the lattice (N^2 points at the end)

  Then in gnuplot:
    plot "out" using 1:2:3:4 with vectors lt 3
*/




/*
 Change fx, fy, gx, gy
*/

//// These conditions admit known analytic solution
//double fx(double x, double y)
//{
//  return -6*x*y*y -2*x*x*x;
//}
//double fy(double x, double y)
//{
//  return 6*x*x*y + 2*y*y*y;
//}
//double gx(double x, double y)
//{
//  return x*x*x*y*y + 1./3.;
//}
//double gy(double x, double y)
//{
//  return -x*x*y*y*y + 2./3.;
//}

// Rotating fluid
double fx(double x, double y)
{
  return 0;
}
double fy(double x, double y)
{
  return 0;
}
double gx(double x, double y)
{
  double res=0;
  if (y==1)//(fabs(y-1.) < 0.01)
  {
    res = 50;
  }
  if (y==0)//(fabs(y) < 0.01)
  {
    res = -50;
  }
  
  return res;
}
double gy(double x, double y)
{
  double res=0;
  if (x==1)//(fabs(x-1.) < 0.1)
  {
    res = -50;
  }
  if (x==0)//(fabs(x) < 0.1)
  {
    res = 50;
  }
  return res;
}


//// Conditions "forming a circle"
//double fx(double x, double y)
//{
//  return 10000;
//}
//double fy(double x, double y)
//{
//  return 0;
//}
//double gx(double x, double y)
//{
//  return 1/2.-x;
//}
//double gy(double x, double y)
//{
//  return 1/2.-y;
//}



/*
  Build the f vector
*/
void constructF(VectorXd & F, unsigned long int N)
{
  for(unsigned long int i=0;i<N*N;++i)
  {
//    std::cerr << i/N +1 << " " << i% N +1 << std::endl; 
    F(i)     = fx((double)(i/N+1)/(double)(N+1), (double)(i%N+1)/(double)(N+1));
    F(i+N*N) = fy((double)(i/N+1)/(double)(N+1), (double)(i%N+1)/(double)(N+1));
  }
  for (unsigned long int i=2*N*N;i<F.rows();++i)
  {
    F(i)=0;
  }
}

/*
  Build boundary conditions with g
*/
void constructBoundaryIndices(VectorXi & bIndices, unsigned long int N)
{
  // bIndices must be 2*(4N+4) dim
  unsigned long int cpt=0;

  // Bottom
  for(unsigned long int i=0;i<N+2;++i)
  {
    bIndices(cpt)=i;
    ++cpt;
  }
  // Left then right alternatively, from bottom to the top
  for(unsigned long int i=1;i<=N;++i)
  {
    bIndices(cpt)=(N+2)*i;
    ++cpt;
    bIndices(cpt)=N+1+(N+2)*i;
    ++cpt;
  }
  // Top
  for(unsigned long int i=0;i<N+2;++i)
  {
    bIndices(cpt)=(N+1)*(N+2)+i;
    ++cpt;
  }
}

/*
  Build the vector containing the indices 
*/
void constructInnerIndices(VectorXi & iIndices, unsigned long int N)
{
  unsigned long int cpt=0;
  for(unsigned long int k=1;k<=N;++k)
  {
    for(unsigned long int j=1;j<=N;++j)
    {
      iIndices(cpt) = k*(N+2)+j;
      ++cpt;
    }
  }
}

//// Analytic solution to the first boundary condition
//void giveSolution(VectorXd & sol, unsigned long int N)
//{
//  double x=0;
//  double y=0;
//  for(unsigned long int i=0;i<N;++i)
//  {
//    x=(double)i/(double)(N+1);
//    for(unsigned long int j=0;j<N;++j)
//    {
//      y=(double)j/(double)(N+1);
//      sol(i*N+j) =  x*x*x*y*y + 1./3.;
//    }
//  }
//  for(unsigned long int i=0;i<N;++i)
//  {
//    x=(double)i/(double)(N+1);
//    for(unsigned long int j=0;j<N;++j)
//    {
//      y=(double)j/(double)(N+1);
//      sol(i*N+j + N*N) = -x*x*y*y*y + 2./3.;
//    }
//  }
////  std::cout << "P constante" << std::endl;
//}











/*      
   Main
   
   
*/

int main(int argc, char ** argv)
{
//  int N=20; // we manipulate unsigned long int indices (from 0 to 4294967295), so N MUST BE < 10920
  int N=20;
  if(argc>1)
    N=atoi(argv[1]);
  
  
  std::cerr << "N " << N << std::endl; 
  
  unsigned long int N2=N*N;
  double h=1./(N+1);
  
  double eps=0.000000001;
  double nu=0.01;
  
  double nu_h2 = nu/(h*h);
  double inv2h = 1./(2.*h);
  
  
  // Boundary indices
  VectorXi bIndices(4*(N+1));
  constructBoundaryIndices(bIndices,N);

  
  // Inner indices
  VectorXi iIndices(N2);
  constructInnerIndices(iIndices,N);
  
  // Build the system A X = f
  
  // A
  SpMat A( 2*N2+(N+1)*(N+1), 2*(N+2)*(N+2)+(N+1)*(N+1));
  unsigned long int nbEntriesA = 2*5*N2 + 2*4*N2 + 2*4*(N+1)*(N+1);

  std::vector<T> listA;
  listA.reserve(nbEntriesA);
  
//  unsigned long int shift=1;
//  unsigned long int shiftP=0;
  
  unsigned long int kk=0; // for the conversion inner indices (0..N*N-1) -> global indices (0..(N+2)*(N+2)-)
  
  // For the velocity, let kk be the global index. Neigbors of the point kk corresponding to the point (i,j):
  // Bottom (i-1,j) : kk-(N+2)
  // Upper  (i+1,j) : kk+(N+2)
  // Left   (i,j-1) : kk-1
  // Right  (i,j+1) : kk+1
  
  // For the pressure, let k be an inner index. Neigbors:
  // Bottom left  (i-1/2,j-1/2) : k + k/N
  // Bottom right (i-1/2,j+1/2) : k + k/N + 1
  // Upper left   (i+1/2,j-1/2) : k + k/N + N+1
  // Upper right  (i+1/2,j+1/2) : k + k/N + N+1 + 1
  
  // Loop on inner indices
  for(unsigned long int k=0;k<iIndices.rows();++k)
  {
    kk = iIndices(k);
    
    // U_x
    listA.push_back(T(k, kk - (N+2)          , -nu_h2  ));
    listA.push_back(T(k, kk + N+2            , -nu_h2  ));
    listA.push_back(T(k, kk - 1              , -nu_h2  ));
    listA.push_back(T(k, kk + 1              , -nu_h2  ));
    listA.push_back(T(k, kk                  , 4*nu_h2  ));
  
    // U_y
    listA.push_back(T(k+N2, kk + (N+2)*(N+2) - (N+2)          , -nu_h2  ));
    listA.push_back(T(k+N2, kk + (N+2)*(N+2) + N+2            , -nu_h2  ));
    listA.push_back(T(k+N2, kk + (N+2)*(N+2) - 1              , -nu_h2  ));
    listA.push_back(T(k+N2, kk + (N+2)*(N+2) + 1              , -nu_h2  ));
    listA.push_back(T(k+N2, kk + (N+2)*(N+2)                  , 4*nu_h2  ));
    
    // P
    listA.push_back(T(k, 2*(N+2)*(N+2) + k + k/N + N+1 + 1 , inv2h  ));
    listA.push_back(T(k, 2*(N+2)*(N+2) + k + k/N + N+1    , inv2h  ));
    listA.push_back(T(k, 2*(N+2)*(N+2) + k + k/N + 1      , -inv2h  ));
    listA.push_back(T(k, 2*(N+2)*(N+2) + k + k/N          , -inv2h  ));
    
    listA.push_back(T(k+N2, 2*(N+2)*(N+2) + k + k/N + N+1 + 1 , inv2h  ));
    listA.push_back(T(k+N2, 2*(N+2)*(N+2) + k + k/N + N+1    , -inv2h  ));
    listA.push_back(T(k+N2, 2*(N+2)*(N+2) + k + k/N + 1      , inv2h  ));
    listA.push_back(T(k+N2, 2*(N+2)*(N+2) + k + k/N          , -inv2h  ));
  }
  

  // For div U = 0:
  
  // Loop on global indices, and we do nothing whan i=N+2 (i.e. when k%(N+2)!=-1), or when j=N+2 (i.e. when k>=(N+2)*(N+1))
  // Let k be a global index. Neigbors:
  // Right       (i,j+1)   : k+1
  // Upper       (i+1,j)   : k+ N+2
  // Upper right (i+1,j+1) : k + N+2 +1

  unsigned long int cpt=0; // will be in 0..(N+1)*(N+1)

  for(unsigned long int k=0;k<(N+2)*(N+1);++k)
  {
    if(k%(N+2) !=(N+1))
    {
      listA.push_back(T(2*N2 + cpt, k           , inv2h  ));
      listA.push_back(T(2*N2 + cpt, k + 1       , inv2h  ));
      listA.push_back(T(2*N2 + cpt, k + N+2     , -inv2h  ));
      listA.push_back(T(2*N2 + cpt, k + N+2 +1  , -inv2h  ));

      listA.push_back(T(2*N2 + cpt, (N+2)*(N+2) + k           , inv2h  ));
      listA.push_back(T(2*N2 + cpt, (N+2)*(N+2) + k + 1       , -inv2h  ));
      listA.push_back(T(2*N2 + cpt, (N+2)*(N+2) + k + N+2     , inv2h  ));
      listA.push_back(T(2*N2 + cpt, (N+2)*(N+2) + k + N+2 +1  , -inv2h  ));
  
      // epsilon part
      listA.push_back(T(2*N2 + cpt, 2*(N+2)*(N+2) + cpt           , -eps  ));

      ++cpt;
    }
  }

  A.setFromTriplets(listA.begin(), listA.end());
//  std::cout << "A:" << std::endl << A << std::endl;
  
  
  
  
  // Build vector f and update it from boundary conditions
  VectorXd F(2*N*N+(N+1)*(N+1));
  constructF(F,N);
//  std::cerr << "F" << std::endl << F << std::endl;
  
  unsigned long int r=0;
  for (unsigned long int k=0;k<bIndices.rows();++k)
  {
    // For each non null element of the kth column of A, multiply it by the corresponding g and put it on the right member
    // U_x
    for (SpMat::InnerIterator it(A,bIndices[k]); it; ++it)
    {
      r = it.row();
      F(r) = F(r) - it.value()*gx((double)(bIndices[k]/(N+2))/(double)(N+1), (double)(bIndices[k]%(N+2))/(double)(N+1));
    }
    
    // U_y
    for (SpMat::InnerIterator it(A,bIndices[k]+(N+2)*(N+2)); it; ++it)
    {
      r = it.row();
      F(r) = F(r) - it.value()*gy((double)(bIndices[k]/(N+2))/(double)(N+1), (double)(bIndices[k]%(N+2))/(double)(N+1));
    }
  }
  

  // Note: DID NOT SUCCEED TO DELETE COLUMNS OF THE SPARSE MATRIX A IN EIGEN. CREATE A NEW MATRIX AA.
  SpMat AA(2*N2+(N+1)*(N+1), 2*N2+(N+1)*(N+1));
  std::vector<T> listAA;
  
  // U_x
  for(unsigned long int k=0;k<iIndices.rows();++k)
  {
    for (SpMat::InnerIterator it(A,iIndices[k]); it; ++it)
    {
      listAA.push_back(T(it.row(), k, it.value() ));
    }
  }
  
  // U_Y
  for(unsigned long int k=0;k<iIndices.rows();++k)
  {
    for (SpMat::InnerIterator it(A,iIndices[k]+(N+2)*(N+2)); it; ++it)
    {
      listAA.push_back(T(it.row(), k+iIndices.rows(), it.value() ));
    }
  }
  
  // P
  for(unsigned long int k=0;k<(N+1)*(N+1);++k)
  {
    for (SpMat::InnerIterator it(A,k+2*(N+2)*(N+2)); it; ++it)
    {
      listAA.push_back(T(it.row(), k+2*iIndices.rows(), it.value() ));
    }
  }
  
  AA.setFromTriplets(listAA.begin(), listAA.end());
//  std::cout << "new A" << std::endl << AA << std::endl;
  
  
  // Solver
  Eigen::SimplicialLDLT<SpMat> SLDLT;
  SLDLT.compute(AA);
  VectorXd X;
  X = SLDLT.solve(F);
  
  // Get the velocity
  VectorXd U(2*N2);
  for(unsigned long int k=0;k<2*N2;++k)
  {
    U(k) = X(k);
  }
  // "normalize"
  U /= (U.norm()*2);
  
  
//  // For the analytic solution of the first example of boundary condition  
//  VectorXd sol(2*N2);
//  giveSolution(sol,N);
//  std::cerr << "L2 Error : " << (U-sol).squaredNorm() << std::endl;
  


  // Display U with the format. Use gnuplot to see the result
  // x y x+u_x y+v_y
  for(unsigned long int k=0;k<N2;++k)
  {
    std::cout << (double)(k/N+1)/(double)(N+1) << " " << (double)(k%N+1)/(double)(N+1) << " " << U(k) << " " << U(k+N2) << std::endl;
  }
}

