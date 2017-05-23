
#ifndef  _BANK_H_    /* only process this file once */
#define  _BANK_H_


#include "Configuration.h"
#include <vector>


struct BankNode{
	std::vector<Configuration> configs;
	BankNode* left;
	BankNode* right;
};

//s


class Bank{

	//This will have to act like a binary tree of configurations, where the ordering is adj matrix in canonical form
	//TODO make this a self-balancing tree
	
public:
	Bank(bool b = false){
		this->isMainBank = b;
		root = new BankNode();
		_size = 0;
	}
	
	~Bank(){
		//This should traverse and delete all the configurations. To be honest,
		// we don't ever actually delete the bank though, so this doesn't need
		// to be implemented
		delete root;
	}
	
	int add(Configuration c);//RETURNS 1 IF ADDED, 0 IF ALREADY IN BANK
	int size();
	void printDetails();
	
	void printDuplicates();
	void recPrintDuplicates(BankNode* node);

	
private:
	bool isMainBank;
	int _size;
	BankNode* root;
	int recursiveAdd(Configuration& c, BankNode* node);
	void recPrint(BankNode* node);

};



#endif
