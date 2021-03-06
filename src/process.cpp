#define MAXN WORDSIZE
#define ONE_WORD_SETS 1
#include "configuration.h"
#include "enumeration.h"
#include <iostream>
#include <fstream>

void readPoints(std::istream& f, Configuration& c) {
  double d;
  int num_of_spheres = c.p.rows() / 3;
  for (int i = 0; i < num_of_spheres; i++) {
    for (int j = 0; j < 3; j++) {
      f >> d;
      c.p(3 * i + j) = d;
    }
  }
  for (int i = 0; i < num_of_spheres - 1; i++) {
    for (int j = i + 1; j < num_of_spheres; j++) {
      double dist = 0;
      for (int k = 0; k < 3; k++) {
        dist += pow(c.p(3 * i + k) - c.p(3 * j + k), 2);
      }
      if (fabs(dist - 1) <= 1e-6) {
        c.addEdge(i, j);
      }
    }
  }
  c.canonize();
}


int main(int argc, char** argv) {
  Enumeration enumeration(true);
  std::cout.precision(16);
  Configuration n[11980];
  Configuration c[11994];
  std::ifstream f;
  f.open("n12.txt");
  if (!f.is_open()) {
    return 1;
  }

  std::ifstream g;
  g.open("parsedoutput.txt");
  if (!g.is_open()) {
    return 1;
  }
  for (int i = 0; i < 11980; i++) {
    readPoints(f, n[i]);
  }
  for (int i = 0; i < 11980; i++) {
    enumeration.add(n[i]);
  }

  for (int i = 0; i < 11994; i++) {
    readPoints(g, c[i]);
  }
  for (int i = 0; i < 11994; i++) {
    enumeration.add(c[i]);
  }

  enumeration.printDuplicates();
  std::cout << enumeration.size() << std::endl;

  return 0;
}
