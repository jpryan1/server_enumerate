#include "configuration.h"
#include <random>

int Configuration::project() {
  VectorXd copy = this->p;
  int temp = project(&copy, &copy);
  if (!temp) return 0;
  this->p = copy;
  return 1;
}


int Configuration::project(VectorXd* old, VectorXd* proj) {
  // TODO include constraint that projection is orthogonal-ish
  // this function takes the vector OLD and solves newtons method on the
  // constraint equations to find a zero.
  VectorXd initial = *old;  // this makes a copy
  double jump_size, F_size;

  ConfigVector delx;

  // Deflation technique - lookup
  // Now allocate
  MatrixXd rigid_x;

  MatrixXd jacob_inv(3 * num_of_spheres, this->num_of_contacts + 6);
  VectorXd F_vec(this->num_of_contacts + 6);
  ConfigVector rand;
  int iterations = 0;

  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0, 1e-12);

  for (int i = 0; i < 4; i++) {
    do {
      iterations++;
      rigid_x = MatrixXd::Zero(this->num_of_contacts + 6,
                               3 * num_of_spheres);  // zero it out
      populateRigidityMatrix(rigid_x,
                             initial);  // repopulate it wrt new initial vector
      rigid_x *= 2;  // now its the jacobian

      populate_F_vec(initial, F_vec);
      F_size = F_vec.norm();
      if (F_size <= NEWTON_TOL) {
        break;
      }

      if (this->isRegular) {
        *proj =  initial - rigid_x.fullPivLu().solve(F_vec);
      } else {
        jacob_inv = pseudoInverse(rigid_x);
        // (rigid_x.transpose()*rigid_x).inverse()*rigid_x;
        // F_vec is constraints evaluated on initial
        *proj = initial - jacob_inv * F_vec;
      }

      delx = *proj - initial;
      jump_size = (delx).norm();
      if (jump_size > DELXMAX) {
        delx.normalize();
        *proj = initial + DELXMAX * delx;
      }

      initial = *proj;
    } while (jump_size > NEWTON_TOL && iterations < MAX_NEWTON_ITERATIONS);

    if (F_size <= NEWTON_TOL) break;

    for (int i = 0; i < 3 * num_of_spheres; i++) {
      rand(i) =  distribution(generator);
    }
    initial = *old + rand;
  }
  populate_F_vec(*proj, F_vec);
  F_size = F_vec.norm();
  // I am not sure why we need to recalculate it, but without
  // recalculation, we were having bad clusters sneak by

  if (F_size > NEWTON_TOL) {
    return 0;
  }
  return 1;
}


void Configuration::populate_F_vec(VectorXd& initial, VectorXd& F_vec) {
  F_vec = MatrixXd::Zero(F_vec.rows(), 1);
  int k = 0;
  for (int i = 0; i < num_of_spheres - 1; i++) {
    for (int j = i + 1; j < num_of_spheres; j++) {
      if (this->hasEdge(i, j)) {
        for (int m = 0; m < 3; m++) {
          F_vec(k) += pow(initial(3 * i + m) - initial(3 * j + m), 2);
        }
        F_vec(k) -= 1;
        k++;
      }
    }
  }
}
