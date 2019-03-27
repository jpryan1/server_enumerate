#define MAXN WORDSIZE
#define ONE_WORD_SETS 1
// The above is especially used since our graphs will have small order,
// see nauty documentation

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <queue>
#include <mutex>
#include <thread>
#include "Configuration.h"
#include "Bank.h"
#define NUM_THREADS 16


std::queue<Configuration> Queue;
std::mutex enumeration_lock;
std::mutex queue_lock;


Configuration initialCluster();
void breakContactsAndAdd(Configuration& current, Enumeration& predim);
void enumerateClusters(Configuration* initial, Enumeration* enumeration);
void debug();

// LATER Consider fixed size vectorization (16 byte alignment) via Eigen
// NOTE: nauty "graph" is just unsigned long (bitwise adj matrix)
// Timer timer;


int main(int argc, char** argv) {
  Eigen::initParallel();
  std::cout.precision(16);

  // Start with a basic cluster formed as a collection of pyramids.
  Configuration c = initialCluster();
  Enumeration* enumeration = new Enumeration(true);

  std::string opts = "\tOptions:\n\t\t(a)nimate\n\t\t(d)ebug\n\t\t(r)un";
  if (argc < 2) {
    std::cout << "usage: ./build/enumerate_clusters {option}\n" << opts <<
              std::endl;
  }
  char option = *(argv[1]);

  if (option == 'a') {
  } else if (option == 'r') {
    enumerateClusters(&c, enumeration);
  } else if (option == 'd') {
    debug();
  }
  return 0;
}


void threadFunc(Enumeration* predim, Enumeration* enumeration) {
  Configuration current;
  do {
    queue_lock.lock();
    current = Configuration(Queue.front());
    Queue.pop();
    queue_lock.unlock();
    int added;
    enumeration_lock.lock();
    if (!current.canonize()) {
      enumeration_lock.unlock();
      continue;
    }
    added = enumeration->add(current);
    enumeration_lock.unlock();
    if (!added) {
      continue;
    }
    breakContactsAndAdd(current, *predim);
  } while (Queue.size() > 0);
}


void enumerateClusters(Configuration* initial, Enumeration* enumeration) {
  double t = time(0);
  Queue.push(*initial);
  Configuration current;
  Enumeration predim;
  while (Queue.size() > 0 && Queue.size() < 50) {
    current = Configuration(Queue.front());
    Queue.pop();
    if (!current.canonize()) continue;
    if (!enumeration->add(current)) {
      // returns 0 if already in enumeration, otherwise adds and returns 1
      continue;
    }
    breakContactsAndAdd(current, predim);
  }
  if (Queue.size() > 0) {
    std::cout << "Calling threads" << std::endl;
    std::thread threads[NUM_THREADS];
    for (int i = 0; i < NUM_THREADS; i++) {
      threads[i] = std::thread(threadFunc, &predim, enumeration);
    }

    for (auto& th : threads) th.join();
  }
  std::cout << "Enumeration is size " << enumeration->size() << std::endl;
  std::cout << "Aux Enumerations are " << predim.size() << std::endl;
  t = time(0) - t;
  std::cout << "Elapsed time: " << t << " seconds." << std::endl;
  enumeration->printDetails();
}


Configuration initialCluster() {
  double points[NUM_OF_SPHERES * 3];

  std::ifstream clusterFile;
  clusterFile.open("first_cluster8.txt");
  if (clusterFile.is_open()) {
    for (int i = 0; i < NUM_OF_SPHERES * 3; i++) {
      clusterFile >> points[i];
    }
    clusterFile.close();
  } else {
    std::cout << "Failed to open initial cluster file!" << std::endl;
    return Configuration();
  }

  graph g[NUM_OF_SPHERES];
  memset(g, 0, NUM_OF_SPHERES * sizeof(graph));

  Configuration* c = new Configuration(points, g);

  for (int i = 0; i < NUM_OF_SPHERES - 3; i++) {
    for (int j = 1; j < 4; j++) {
      c->addEdge(i, i + j);
    }
  }
  c->addEdge(NUM_OF_SPHERES - 3, NUM_OF_SPHERES - 2);
  c->addEdge(NUM_OF_SPHERES - 3, NUM_OF_SPHERES - 1);
  c->addEdge(NUM_OF_SPHERES - 2, NUM_OF_SPHERES - 1);

  c->canonize();
  return *c;
}


void breakContactsAndAdd(Configuration& current, Enumeration& predim) {
  int dim;
  Configuration copy;
  std::vector<Configuration> walkedTo;

  // These loops iterate through all edges in the graph of current
  for (int i = 0; i < NUM_OF_SPHERES; i++) {
    for (int j = i + 1; j < NUM_OF_SPHERES; j++) {
      if (!current.hasEdge(i, j)) {
        continue;
      }
      // We make a copy, delete an edge of that copy, then either enter that
      // copy into the queue, or delete it, depending on whether it is rigid.
      copy = Configuration(current);
      copy.deleteEdge(i, j);
      copy.ptype.perms.clear();
      enumeration_lock.lock();
      if (!copy.canonize()) {
        enumeration_lock.unlock();
        continue;
      }
      if (!predim.add(copy)) {
        enumeration_lock.unlock();
        continue;
      }

      enumeration_lock.unlock();
      dim = copy.dimensionOfTangentSpace(true);


      if (dim == 0) {
        breakContactsAndAdd(copy, predim);
      } else if (dim == 1) {
        walkedTo = copy.walk();

        for (int k = 0; k < walkedTo.size(); k++) {  // walking can fail!
          int temp = walkedTo[k].dimensionOfTangentSpace(false);
          if (temp > 0) continue;

          queue_lock.lock();
          Queue.push(walkedTo[k]);
          queue_lock.unlock();
        }
      }
    }
  }
}


void debug() {
  Configuration first;
  std::ifstream debugFile;
  debugFile.open("debug.txt");
  if (debugFile.is_open()) {
    first.readClusterFromFile(debugFile);

    debugFile.close();
  } else {
    std::cout << "Failed to open debug file!" << std::endl;
  }
}
