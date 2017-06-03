
#include "Bank.h"

int Bank::add(Configuration c){ //RETURNS 1 IF ADDED, 0 IF ALREADY IN BANK

	if(root->configs.size()<1){
		_size++;
		root->configs.push_back(c);
		return 1;
	}
	int temp = recursiveAdd(c, root);
	if(temp){
		_size++;
	}return temp;
}

int Bank::size(){
	return _size;
}
int Bank::recursiveAdd(Configuration& c, BankNode* node){

	int comp = c.compareGraph(node->configs[0]);
	if(comp == 0){
		
		for(int i=0; i<node->configs.size(); i++){
			if(c.permMatches(node->configs[i])){
				return 0;
			}
		}
	/*	if(isMainBank){
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
	}else if(comp < 0){
		if(!node->left){
			node->left = new BankNode();
			node->left->configs.push_back(c);
			return 1;
		}else{
			return recursiveAdd(c, node->left);
		}
	}else{
		if(!node->right){
			node->right = new BankNode();
			node->right->configs.push_back(c);
			return 1;
		}else{
			return recursiveAdd(c, node->right);
		}
		
	}
	
}



void Bank::printDuplicates(){
	if(!root){
		return;
	}
	if(root->configs.size()<1){
		return;
	}
	recPrintDuplicates(root);
}
void Bank::recPrintDuplicates(BankNode* node){
	
	if(node->configs.size()<=1){}
	else{
		for(int i=0; i<node->configs.size(); i++) node->configs[i].printDetails();
		std::cout<<"\n\n\n"<<std::endl;
	}
	if(!node->left){}
	else{
		recPrintDuplicates(node->left);
	}
	if(!node->right){}
	else{
		recPrintDuplicates(node->right);
	}
}








void Bank::printDetails(){
	//std::cout<<"Printing bank details..."<<std::endl;
	if(!root){
		//std::cout<<"No root!"<<std::endl;
		return;
	}
	if(root->configs.size()<1){
		//std::cout<<"Empty root!"<<std::endl;
		return;
	}
	//std::cout<<"Starting with root..."<<std::endl;
	recPrint(root);
	//std::cout<<"Done printing bank details!"<<std::endl;
	
}
void Bank::recPrint(BankNode* node){
	
	if(node->configs.size()<1){
		//std::cout<<"This node has no configs!"<<std::endl;
	}else{
		for(int i=0; i<node->configs.size(); i++) node->configs[i].printDetails();
		//std::cout<<"This node has "<<node->configs.size()<<" configs"<<std::endl;
	}
	if(!node->left){
		//std::cout<<"No left child!"<<std::endl;
	}else{
		//std::cout<<"Going left..."<<std::endl;
		recPrint(node->left);
	}
	if(!node->right){
		//std::cout<<"No right child!"<<std::endl;
	}else{
		//std::cout<<"Going right..."<<std::endl;
		recPrint(node->right);
	}
	//std::cout<<"Going back up..."<<std::endl;
}



