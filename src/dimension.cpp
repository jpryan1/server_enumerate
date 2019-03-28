#include "configuration.h"

MatrixXd Configuration::getRightNullSpace(MatrixXd rigid, bool* null_flag) {
  if (this->isRegular) {
    FullPivLU<MatrixXd> rightlu(rigid);
    MatrixXd right_null_space = rightlu.kernel();
    if (right_null_space.cols() == 1
        && right_null_space.isZero(1e-5)) {
      // 1e-5 is the precision with which to check.
      *null_flag = true;
    }
    return right_null_space;
  }
  JacobiSVD<MatrixXd> svd(rigid, ComputeFullU | ComputeFullV);
  MatrixXd svdv = svd.matrixV();
  MatrixXd sing = svd.singularValues();
  int smallsing = 0;
  for (int i = 0; i < sing.size(); i++) {
    if (sing(i) < 1e-6) {
      smallsing = sing.size() - i;
    }
  }
  int V = smallsing + svdv.cols() - sing.size();
  if (V == 0) {
    *null_flag = true;
  }
  return svdv.rightCols(V);
}


int Configuration::dimensionOfTangentSpace(bool useNumericalMethod) {
  // Construct rigidity matrix
  // Rigidity matrix has dimension numContacts x 3*numofspheres
  // First we calculate dimensions of matrix
  // - rigidity matrix has dimension numContacts x 3*numofspheres

  // Allocate
  MatrixXd rigid_x = MatrixXd::Zero(this->num_of_contacts + 6,
                                    3 * num_of_spheres);

  // Populate
  populateRigidityMatrix(rigid_x, this->p);

  this->isRegular = false;

  // Find right null space
  bool null_flag = false;
  MatrixXd right_null_space = getRightNullSpace(rigid_x, &null_flag);
  if (null_flag) {
    if (0 <= 3 * num_of_spheres - 6 - num_of_contacts) this->isRegular = true;
    return 0;
  }
  int V = right_null_space.cols();
  this->v = right_null_space.col(0);
  this->v = this->v / this->v.norm();
  if (V <= 3 * num_of_spheres - 6 - num_of_contacts) {
    this->isRegular = true;
    return V;
  } else {
    this->isRegular = false;
  }
  // Find left null space
  null_flag = false;
  MatrixXd left_null_space = getRightNullSpace(rigid_x.transpose(), &null_flag);
  int W;
  if (null_flag) {
    W = 0;
  } else {
    W = left_null_space.cols();
  }
  if (W == 0) {
    if (useNumericalMethod) {
      return numerical_findDimension(right_null_space);
    }
    return V;
  }

  // Compute Q matrices, check for sign-definiteness via eigendecomposition
  MatrixXd Q = MatrixXd::Zero(V, V);
  MatrixXd R_vi = MatrixXd::Zero(this->num_of_contacts + 6, 3 * num_of_spheres);
  VectorXd vi(3 * num_of_spheres);
  VectorXd eigs(V);
  bool flag;
  // this can be vectorized better... TODO
  for (int k = 0; k < W; k++) {
    for (int i = 0; i < V; i++) {
      vi << right_null_space.col(i);
      rigid_x = MatrixXd::Zero(this->num_of_contacts + 6, 3 * num_of_spheres);
      populateRigidityMatrix(R_vi, vi);
      for (int j = 0; j < V; j++) {
        Q(i, j) = left_null_space.col(k).transpose() *
                  R_vi * right_null_space.col(j);
      }

      // test Q for sign definiteness here.
      SelfAdjointEigenSolver<MatrixXd> es(Q, EigenvaluesOnly);
      eigs = es.eigenvalues();
      flag = true;
      if (eigs(0) > 0) {
        for (int idx = 0; idx < V; idx++) {
          if (eigs(idx) <= 1e-5) {
            flag = false;
          }
        }
        if (flag) {
          return 0;
        }
      } else if (eigs(0) < 0) {
        for (int idx = 0; idx < V; idx++) {
          if (eigs(idx) >= -1e-5) {
            flag = false;
          }
        }
        if (flag) {
          return 0;
        }
      }
      vi = Matrix<double, 3 * NUM_OF_SPHERES, 1>::Zero();
      // zero out for next iteration
      R_vi = MatrixXd::Zero(this->num_of_contacts + 6, 3 * num_of_spheres);
    }
  }
  if (useNumericalMethod) {
    return numerical_findDimension(right_null_space);
  }
  return V;
}


int Configuration::numerical_findDimension(MatrixXd& right_null_space) {
  std::vector<VectorXd> basis;
  VectorXd jump(3 * num_of_spheres), proj(3 * num_of_spheres),
    tang(3 * num_of_spheres), orth(3 * num_of_spheres);
  for (int i = 0; i < right_null_space.cols(); i++) {
    jump = DEL_S0 * right_null_space.col(i) + this->p;

    if (project(&jump, &proj)) {
      if ((proj - jump).norm() <= TOLMAX) {
        tang = proj - this->p;
        if (tang.norm() < TOLMIN) {
          continue;
        }
        orth = tang;
        for (int j = 0; j < basis.size(); j++) {
          orth = orth - orth.dot(basis[j]) * basis[j];
        }
        if (orth.norm() > vTol) {
          basis.push_back(orth / orth.norm());
        }
      }
    }
    jump = -DEL_S0 * right_null_space.col(i) + this->p;
    if (project(&jump, &proj)) {
      if ((proj - jump).norm() <= TOLMAX) {
        tang = proj - this->p;
        if (tang.norm() < TOLMIN) {
          continue;
        }
        orth = tang;
        for (int j = 0; j < basis.size(); j++) {
          orth = orth - orth.dot(basis[j]) * basis[j];
        }
        if (orth.norm() > vTol) {
          basis.push_back(orth / orth.norm());
        }
      }
    }
  }
  if (basis.size() > 0) {
    this->v = basis[0];
    this->v = this->v / this->v.norm();
  }
  return basis.size();
}
