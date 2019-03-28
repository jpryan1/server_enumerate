#include "enumeration.h"

int Enumeration::add(Configuration c) {
  if (root->configs.size() < 1) {
    _size++;
    root->configs.push_back(c);
    return 1;
  }
  int temp = recursiveAdd(c, root);
  if (temp) {
    _size++;
  } return temp;
}


int Enumeration::size() {
  return _size;
}


int Enumeration::recursiveAdd(Configuration& c, EnumerationNode* node) {
  int comp = c.compareGraph(node->configs[0]);
  if (comp == 0) {
    for (int i = 0; i < node->configs.size(); i++) {
      if (c.permMatches(node->configs[i])) {
        return 0;
      }
    }
    /*	if(isMainEnumeration){
    		std::cout<<"Compared the following, no match"<<std::endl;
    		c.printDetails();
    		node->configs[0].printDetails();
    		c.permMatches(node->configs[0],true);
    //			MatrixXd F_vec(c.num_of_contacts+6, 1);
    //			c.populate_F_vec(c.p, F_vec);
    //			MatrixXd F_vec2(node->configs[0].num_of_contacts+6, 1);
    //			node->configs[0].populate_F_vec(node->configs[0].p, F_vec2);
    //			std::cout<<"Fs "<<F_vec.norm()<<" "<<F_vec2.norm();
    //			std::cout<<"Orb distance:"<<std::endl;
    //			c.orbitMatches(node->configs[0],0, true);
    //			std::cout<<"Orbs"<<std::endl;
    //			c.printOrbits();
    //			node->configs[0].printOrbits();
    //			std::cout<<"\n\n\n\n"<<std::endl;
    	}*/
    node->configs.push_back(c);
    return 1;
  } else if (comp < 0) {
    if (!node->left) {
      node->left = new EnumerationNode();
      node->left->configs.push_back(c);
      return 1;
    } else {
      return recursiveAdd(c, node->left);
    }
  } else {
    if (!node->right) {
      node->right = new EnumerationNode();
      node->right->configs.push_back(c);
      return 1;
    } else {
      return recursiveAdd(c, node->right);
    }
  }
}


void Enumeration::printDuplicates() {
  if (!root) {
    return;
  }
  if (root->configs.size() < 1) {
    return;
  }
  recPrintDuplicates(root);
}


void Enumeration::recPrintDuplicates(EnumerationNode* node) {
  if (node->configs.size() <= 1) {
    // TODO no-op?
  } else {
    for (int i = 0; i < node->configs.size(); i++) {
      node->configs[i].printDetails();
    }
    std::cout << "\n\n\n" << std::endl;
  }
  if (!node->left) {
    // TODO no-op?
  } else {
    recPrintDuplicates(node->left);
  }
  if (!node->right) {
    // TODO no-op?
  } else {
    recPrintDuplicates(node->right);
  }
}


void Enumeration::printDetails() {
  if (!root) {
    return;
  }
  if (root->configs.size() < 1) {
    return;
  }
  recPrint(root);
}


void Enumeration::recPrint(EnumerationNode* node) {
  if (node->configs.size() < 1) {
  } else {
    for (int i = 0; i < node->configs.size(); i++) {
      node->configs[i].printDetails();
    }
  }
  if (!node->left) {
  } else {
    recPrint(node->left);
  }
  if (!node->right) {
  } else {
    recPrint(node->right);
  }
}
