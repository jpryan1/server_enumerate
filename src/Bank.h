#ifndef  _ENUMERATION_H_    /* only process this file once */
#define  _ENUMERATION_H_

#include "Configuration.h"
#include <vector>


struct EnumerationNode {
  std::vector<Configuration> configs;
  EnumerationNode* left;
  EnumerationNode* right;
};


class Enumeration {
  // This class will act like a binary tree of configurations, where the
  // ordering is adj matrix in canonical form
 public:
  Enumeration(bool b = false) {
    this->isMainEnumeration = b;
    root = new EnumerationNode();
    _size = 0;
  }

  ~Enumeration() {
    // Should traverse and delete all the configurations. TODO implement
    // this.
    delete root;
  }

  int add(Configuration c);  // Returns 1 if added, 0 if already in enumeration
  int size();
  void printDetails();

  void printDuplicates();
  void recPrintDuplicates(EnumerationNode* node);

 private:
  bool isMainEnumeration;
  int _size;
  EnumerationNode* root;
  int recursiveAdd(Configuration& c, EnumerationNode* node);
  void recPrint(EnumerationNode* node);
};



#endif
