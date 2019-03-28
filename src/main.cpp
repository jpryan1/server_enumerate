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
#include "configuration.h"
#include "enumeration.h"
#define NUM_THREADS 16


std::queue<Configuration> Queue;
std::mutex enumeration_lock;
std::mutex queue_lock;


// LATER Consider fixed size vectorization (16 byte alignment) via Eigen
// NOTE: nauty "graph" is just unsigned long (bitwise adj matrix)


Configuration initialCluster(int num_of_spheres) {
  double points[num_of_spheres * 3];
  std::ifstream clusterFile;
  clusterFile.open("first_cluster.txt");
  if (clusterFile.is_open()) {
    for (int i = 0; i < num_of_spheres * 3; i++) {
      clusterFile >> points[i];
    }
    clusterFile.close();
  } else {
    std::cout << "Failed to open initial cluster file!" << std::endl;
    return Configuration();
  }

  graph g[num_of_spheres];
  memset(g, 0, num_of_spheres * sizeof(graph));
  Configuration* c = new Configuration(points, g, num_of_spheres);

  for (int i = 0; i < num_of_spheres - 3; i++) {
    for (int j = 1; j < 4; j++) {
      c->addEdge(i, i + j);
    }
  }
  c->addEdge(num_of_spheres - 3, num_of_spheres - 2);
  c->addEdge(num_of_spheres - 3, num_of_spheres - 1);
  c->addEdge(num_of_spheres - 2, num_of_spheres - 1);

  c->canonize();
  return *c;
}


void break_contacts_and_add(Configuration& current, Enumeration& predim) {
  int dim;
  int num_of_spheres = current.p.rows() / 3;
  Configuration copy;
  std::vector<Configuration> walkedTo;
  
  // These loops iterate through all edges in the graph of current
  for (int i = 0; i < num_of_spheres; i++) {
    for (int j = i + 1; j < num_of_spheres; j++) {
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
        break_contacts_and_add(copy, predim);
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


void process_queue(Enumeration* predim, Enumeration* enumeration) {
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
    break_contacts_and_add(current, *predim);
  } while (Queue.size() > 0);
}


int main(int argc, char** argv) {
  Eigen::initParallel();
  std::cout.precision(16);

  int num_of_spheres = std::atoi(argv[1]);

  Configuration initial = initialCluster(num_of_spheres);

  Enumeration* enumeration = new Enumeration(true);

  double t = time(0);
  Queue.push(initial);
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
    break_contacts_and_add(current, predim);
  }
  if (Queue.size() > 0) {
    std::cout << "Calling threads" << std::endl;
    std::thread threads[NUM_THREADS];
    for (int i = 0; i < NUM_THREADS; i++) {
      threads[i] = std::thread(process_queue, &predim, enumeration);
    }

    for (auto& th : threads) th.join();
  }
  std::cout << "Enumeration is size " << enumeration->size() << std::endl;
  std::cout << "Aux Enumerations are " << predim.size() << std::endl;
  t = time(0) - t;
  std::cout << "Elapsed time: " << t << " seconds." << std::endl;
  // enumeration->printDetails();
  return 0;
}
