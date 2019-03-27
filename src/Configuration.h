#ifndef  _CONFIGURATION_H_    /* only process this file once */
#define  _CONFIGURATION_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <ctime>
#include <limits>
#include <algorithm>
#include "Eigen/Dense"
#include "nauty.h"
#include "naugroup.h"

#define NUM_OF_SPHERES 8

#define MAXP 1e5   /* maximum number of permutations to count */

#define TOLMAX 10*DEL_S
#define TOLMIN DEL_S/8
#define vTol 2*DEL_S
#define DEL_S0 5e-2
#define DEL_S 5e-3
#define NEWTON_TOL 8e-16
#define tolA 1e-3
// new tolerance for jump back 1e-3
#define tolD 1e-5
#define DELXMAX 0.02
// tolD should be like 1e-6 ... this may be because we don't baby step yet
#define MAX_NEWTON_ITERATIONS 400

using namespace Eigen;

typedef Matrix<double, 3 * NUM_OF_SPHERES, 1> ConfigVector;

typedef std::pair <int, int> Contact;

typedef struct {
  int np;
  std::vector<int> perms;
} permtype;


template<typename _Matrix_Type_>
// This is taken from https://fuyunfei1.gitbooks.io/c-tips/content/pinv_with_eigen.html
_Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a,
                            double epsilon =
                              std::numeric_limits<double>::epsilon()) {
  Eigen::JacobiSVD<_Matrix_Type_> svd(a, Eigen::ComputeThinU
                                      | Eigen::ComputeThinV);
  double tolerance = epsilon * std::max(a.cols(), a.rows())
                     * svd.singularValues().array().abs()(0);
  return svd.matrixV() * (svd.singularValues().array().abs()
                          > tolerance).select(
           svd.singularValues().array().inverse(), 0).matrix().asDiagonal()
         * svd.matrixU().adjoint();
}


class Configuration {
 public:
  Configuration() {
    this->num_of_contacts = 0;
    memset(this->g, 0, sizeof(graph)*NUM_OF_SPHERES);
  }
  Configuration(double* points, graph* adj);
  ~Configuration() {}

  Configuration(const Configuration& other) {
    num_of_contacts = other.num_of_contacts;
    memcpy(&p, &other.p, 3 * NUM_OF_SPHERES * sizeof(double));
    memcpy(g, other.g, NUM_OF_SPHERES * sizeof(graph));
    memcpy(triangle, other.triangle, 3 * sizeof(int));
    ptype.np = other.ptype.np;
    ptype.perms = other.ptype.perms;
  }


  void chooseTriangle();
  int fixTriangle();
  int checkTriangle();
  void printDetails();


  int compareGraph(Configuration& other);  // returns true if graphs match

  int matchesHelper(Configuration& other,  bool det = false);
  int permMatches(Configuration& other,  bool det = false);

  int canonize();

  void deleteEdge(int i, int j);
  void addEdge(int i, int j);

  int hasEdge(int i, int j);

  int dimensionOfTangentSpace(bool useNumericalMethod = true);
  int numerical_findDimension(MatrixXd& right_null_space);
  MatrixXd getRightNullSpace(MatrixXd rigid, bool* null_flag);

  int project(ConfigVector& old, ConfigVector& proj);
  int project();
  void populate_F_vec(ConfigVector& initial, MatrixXd& F_vec);

  std::vector<Contact> checkForNewContacts(ConfigVector proj,
      bool smallTol = false);
  void populateRigidityMatrix(MatrixXd& rigid, ConfigVector& points);
  std::vector<Configuration> walk();

  ConfigVector getP();
  graph* getG();
  void setTriangle(int* t) {
    memcpy(triangle, t, 3 * sizeof(int));
  }
  int num_of_contacts;
  void printTriangle();
  void printAdj();


  void show(int a);

  void readClusterFromFile(std::istream& file);

  void swap(int i, int k);

  ConfigVector p;
  permtype ptype;
  int triangle[3];

 private:
  bool isRegular;
  ConfigVector v;

  // contains 3*n doubles for points in space

  graph g[NUM_OF_SPHERES];
};

#endif
