#include "configuration.h"


Configuration::Configuration(double* points, graph* adj, int n) {
  
  p = VectorXd(3*n);
  num_of_spheres = n;
  for(int i=0; i<num_of_spheres * 3; i++){
    p(i) = points[i];
  }
  memcpy(g, adj, num_of_spheres * sizeof(graph));
  this->num_of_contacts = 0;
  for (int i = 0; i < num_of_spheres - 1; i++) {
    for (int j = i + 1; j < num_of_spheres; j++) {
      if (this->hasEdge(i, j)) {
        this->num_of_contacts++;
      }
    }
  }
}


void Configuration::printTriangle() {
  for (int i = 0; i < 3; i++) {
    std::cout << triangle[i] << ": ";
    for (int j = 0; j < 3; j++) {
      std::cout << p(3 * triangle[i] + j) << " ";
    } std::cout << std::endl;
  }



  double dist = 0;
  for (int i = 0; i < 3; i++) {
    dist += pow(p(3 * triangle[0] + i) - p(3 * triangle[1] + i), 2);
  } std::cout << "0 1: " << dist << std::endl;
  dist = 0;

  for (int i = 0; i < 3; i++) {
    dist += pow(p(3 * triangle[0] + i) - p(3 * triangle[2] + i), 2);
  } std::cout << "0 2: " << dist << std::endl;
  dist = 0;

  for (int i = 0; i < 3; i++) {
    dist += pow(p(3 * triangle[1] + i) - p(3 * triangle[2] + i), 2);
  } std::cout << "1 2: " << dist << std::endl;
  dist = 0;
}


void Configuration::show(int a) {
  double t;
  t = time(0);
  while (time(0) - t < a) {}
}


void Configuration::deleteEdge(int a, int b) {
  DELELEMENT(g + a, b);
  DELELEMENT(g + b, a);
  this->num_of_contacts--;
}


void Configuration::addEdge(int a, int b) {
  ADDONEEDGE(g, a, b, 1);
  this->num_of_contacts++;
}


int Configuration::hasEdge(int i, int j) {
  return ((g[i] & bit[j]) > 0 || (g[j] & bit[i]) > 0);
}


// TODO the following doesn't account for whether the fixed points are collinear
void Configuration::populateRigidityMatrix(MatrixXd& rigid, VectorXd& x) {
  int row = 0;
  for (int i = 0; i < num_of_spheres - 1; i++) {
    for (int j = i + 1; j < num_of_spheres; j++) {
      if (this->hasEdge(i, j)) {
        // ith trio should be p_i - p_j, jth trio should be p_j-p_i
        for (int k = 0; k < 3; k++) {
          rigid(row, 3 * i + k) = x(3 * i + k) - x(3 * j + k);
          rigid(row, 3 * j + k) = x(3 * j + k) - x(3 * i + k);
        }
        row++;
      }
    }
  }
  rigid(row++, 3 * triangle[0]) = 1;
  rigid(row++, 3 * triangle[0] + 1) = 1;
  rigid(row++, 3 * triangle[0] + 2) = 1;

  rigid(row++, 3 * triangle[1] + 1) = 1;
  rigid(row++, 3 * triangle[1] + 2) = 1;
  rigid(row++, 3 * triangle[2] + 2) = 1;
}


std::vector<Contact> Configuration::checkForNewContacts(ConfigVector proj,
    bool smallTol) {
  std::vector<Contact> newContacts;
  double dist, tol;
  if (smallTol) {
    tol = tolA * 0.001;
  } else {
    tol = tolA;
  }
  for (int i = 0; i < num_of_spheres - 1; i++) {
    for (int j = i + 1; j < num_of_spheres; j++) {
      if (this->hasEdge(i, j)) {
        continue;
      }

      dist = 0;
      for (int k = 0; k < 3; k++) {
        dist += pow(proj(3 * i + k) - proj(3 * j + k) , 2);
      }
      if (dist - 1 <= tol) {
        newContacts.push_back(Contact(i, j));
      }
    }
  }
  return newContacts;
}


void Configuration::printDetails() {
  for (int i = 0; i < num_of_spheres * 3; i += 3) {
    for (int j = 0; j < 3; j++) {
      std::cout << (this->p)(i + j) << " ";
    }
  }
  std::cout << std::endl;
}


void Configuration::readClusterFromFile(std::istream& file) {
  memset(this->g, 0, num_of_spheres * sizeof(graph));
  double d;
  for (int i = 0; i < num_of_spheres; i++) {
    for (int j = 0; j < 3; j++) {
      file >> d;
      this->p(3 * i + j) = d;
    }
  }
  for (int i = 0; i < num_of_spheres - 1; i++) {
    for (int j = i + 1; j < num_of_spheres; j++) {
      double dist = 0;
      for (int k = 0; k < 3; k++) {
        dist += pow(this->p(3 * i + k) - this->p(3 * j + k), 2);
      }
      if (fabs(dist - 1) <= 1e-6) {
        this->addEdge(i, j);
      }
    }
  }

  this->canonize();
}


ConfigVector Configuration::getP() {
  return this->p;
}


graph* Configuration::getG() {
  return this->g;
}
