#include "Configuration.h"
#define MAXP 1e5   /* maximum number of permutations to count */







void countperms(int *p, int n, int *abort, void *userptr)  {
	permtype *ptype;
	ptype = (permtype*)userptr;
	int i;
	if(ptype->np >= MAXP) {
		*abort = 1;
	}
	else {
		ptype->np = ptype->np + 1;
		for(i=0; i < NUM_OF_SPHERES; i++){
			ptype->perms.push_back(p[i]);
		}
	}
}


int Configuration::compareGraph(Configuration& other){
	if(this->num_of_contacts!=other.num_of_contacts){
		if(this->num_of_contacts>other.num_of_contacts) return 1;
		return -1;
	}
	
	setword* g1 = (setword*) this->g;
	setword* g2 = (setword*) other.g;
	
	
	for(int i=0; i<NUM_OF_SPHERES; i++){
		if(g1[i]<g2[i]){
			return -1;
		}else if(g1[i]>g2[i]){
			return 1;
		}
	}return 0;
}

int Configuration::matchesHelper(Configuration& other, bool det){
	if(!fixTriangle()) return 1;
	ConfigVector copy = other.p;
	ConfigVector diff = copy - this->p;
	double temp = diff.norm();
	if(det) std::cout<<temp<<std::endl;
	if(diff.norm()<tolD){
		//this->p = pre_fix;
		return 1;
	}
	for(int i=2; i <3*NUM_OF_SPHERES; i+=3) copy(i)  = -copy(i);
	diff = copy - this->p;
	//this->p = pre_fix;
	if(det) std::cout<<diff.norm()<<std::endl;
	if(diff.norm()<tolD) return 1;
	return 0;
	
}


int Configuration::permMatches(Configuration& other,  bool det){
	ConfigVector holder = this->p;
	for(int i=0; i<ptype.np; i++){
		
		
		for(int j=0; j<NUM_OF_SPHERES; j++){
			int perm = ptype.perms[i*NUM_OF_SPHERES+j];
			for(int k=0; k<3; k++){
				this->p(3*j+k) = holder(3*perm + k);
			}
		}
		//check triangle, if its bad then move on
		if(!checkTriangle()){
			this->p = holder;
			continue;
		}
		bool flag = false;
		if(matchesHelper(other, det)){
			flag = true;
		}
		this->p = holder;
		
		if(flag){
			return 1;
		}
	}
	return 0;
}



int Configuration::checkTriangle(){
	double dist = 0;
	for(int i=0; i<3; i++){
		dist+= pow(p(3*triangle[0]+i)-p(3*triangle[1]+i),2);
	}
	if(fabs(dist-1)>1e-5){
		return 0;
	}
	dist = 0;
	
	for(int i=0; i<3; i++){
		dist+= pow(p(3*triangle[0]+i)-p(3*triangle[2]+i),2);
	}
	if(fabs(dist-1)>1e-5){
		return 0;
	}
	dist = 0;
	
	for(int i=0; i<3; i++){
		dist+= pow(p(3*triangle[2]+i)-p(3*triangle[1]+i),2);
	}
	if(fabs(dist-1)>1e-5){
		return 0;
	}
	return 1;
	
}

void Configuration::chooseTriangle(){
	
	std::vector<int> degreeList[NUM_OF_SPHERES];
	
	for(int i=0; i<NUM_OF_SPHERES; i++){
		int numEdges = 0;
		graph temp = this->g[i];
		while(temp>0){
			if(temp%2==1){
				numEdges++;
			}
			temp = temp>>1;
		}
		degreeList[numEdges].push_back(i);
	}
	
	
	
	//so now the ith slot in degreeList is a vector containing the vertices with degree i
	
	int SBD[NUM_OF_SPHERES];
	//short for SORTED BY DEGREE
	
	int idx = NUM_OF_SPHERES-1;
	for(int i=NUM_OF_SPHERES-1; i>=0; i--){
		for(int j = 0; j<degreeList[i].size(); j++){
			SBD[idx--] = degreeList[i][j];
		}
	}
	
	//now sorted by degree has the vertices sorted by degree
	
	
	//note that the below forloop could take NUM_OF_SPHERES^3 time, but
	//should usually take much less. We are looking among high degree vertices
	//for a triangle.
	for(int i=NUM_OF_SPHERES-1; i>=0; i--){
		for(int j = i-1; j>=0; j--){
			if(this->hasEdge(SBD[i],SBD[j])){
				for(int k=j-1; k>=0; k--){
					if(this->hasEdge(SBD[i],SBD[k])&&this->hasEdge(SBD[j],SBD[k])){
						triangle[0] = SBD[i];
						triangle[1] = SBD[j];
						triangle[2] = SBD[k];
						return;
					}
				}
			}
		}
	}
	std::cout<<"Error: no triangle found."<<std::endl;
	exit(0);
}


int Configuration::fixTriangle(){

	double translate[3];
	for(int i=0; i<3; i++){
		translate[i] = -this->p(3*triangle[0]+i);
	}
	for(int i=0; i<NUM_OF_SPHERES; i++){
		for(int j=0; j<3; j++){
			this->p(3*i+j) += translate[j];
		}
	}
	
	double spherepos1[3];
	for(int i=0; i<3; i++){
		spherepos1[i] = this->p(3*triangle[1]+i);
	}
	Vector3d rotate1_axis;
	//take dot with x axis
	//cos for angle
	
	Vector3d temp;
	
	
	double rotate1_ang;
	
	
	if(fabs(-1-spherepos1[0])<1e-14){ //if we are close to (-1,0,0)
		rotate1_ang = M_PI;
		rotate1_axis << 0,1,0;
		AngleAxisd rotationMatrix1( rotate1_ang, rotate1_axis);
		for(int i=0; i<NUM_OF_SPHERES; i++){
			for(int j=0; j<3; j++){
				temp(j) = this->p(3*i+j);
			}
			temp = rotationMatrix1*temp;
			for(int j=0; j<3; j++){
				this->p(3*i+j) = temp(j);
			}
		}

		
	}
	else if(fabs(1-spherepos1[0])>=1e-14){
		rotate1_ang = acos(spherepos1[0]);
		//cross with x axis for rot axis, normalize
		rotate1_axis << 0, spherepos1[2], -spherepos1[1];
		rotate1_axis.normalize();
		//rotate by negative that angle
		AngleAxisd rotationMatrix1( rotate1_ang, rotate1_axis);
		for(int i=0; i<NUM_OF_SPHERES; i++){
			for(int j=0; j<3; j++){
				temp(j) = this->p(3*i+j);
			}
			temp = rotationMatrix1*temp;
			for(int j=0; j<3; j++){
				this->p(3*i+j) = temp(j);
			}
		}
	}
	
	
	double rotate2_ang;
	
	double spherepos2[3];
	
	for(int i=0; i<3; i++){
		spherepos2[i] = this->p(3*triangle[2]+i);
	}
	
	double dot_temp =spherepos2[1]/(sqrt(   pow(spherepos2[1],2) + pow(spherepos2[2],2)     ));
	
	if(fabs(-1-dot_temp)<1e-14){
		rotate2_ang = M_PI;
	}
	else{
		rotate2_ang = acos(dot_temp);
	}
	if(spherepos2[2]<0){
		rotate2_ang *=-1;
	}
	if(fabs(rotate2_ang)>1e-14){
		AngleAxisd rotationMatrix2;
		rotationMatrix2 = AngleAxisd(-rotate2_ang,  Vector3d::UnitX());
		for(int i=0; i<NUM_OF_SPHERES; i++){
			for(int j=0; j<3; j++){
				temp(j) = this->p(3*i+j);
			}
			temp = rotationMatrix2*temp;
			for(int j=0; j<3; j++){
				this->p(3*i+j) = temp(j);
			}
		}
	}
//	
//	if(this->p(3*triangle[2])<0){
//		std::cout<<" Fixing just failed, what a horrible function"<<std::endl;
//		std::cout<<"Previously: "<<std::endl;
//		cop.printTriangle();
//		std::cout<<"Now"<<std::endl;
//		printTriangle();
//	}
	
	
	
	
	return project();
}


int Configuration::canonize(){
	
	
	//These arrays are passed to sparsenauty in canonization, so that function can write its data somewhere
	//TODO edit nauty to not record data into these
	//long story short - nauty does more things than just canonize, which we don't care about
	//sure, this allocation is done every time this function is called,
	//but I think the work is done at compiletime, not runtime
	//besides, the alternative is passing pointers to arrays to this method, which is ugly.
	int lab[NUM_OF_SPHERES];
	int ptn[NUM_OF_SPHERES];
	int orbits[NUM_OF_SPHERES];
	grouprec *group;
	statsblk stats;
	graph canonized[NUM_OF_SPHERES];

	DEFAULTOPTIONS_GRAPH(options);
	options.getcanon = true;
	
	
	
	/* The following case nauty to call two procedures which store the
	 * group information as nauty runs */
	options.userautomproc = groupautomproc;
	options.userlevelproc = grouplevelproc;
	
	
	

	densenauty(this->g, lab, ptn, orbits, &options, &stats, 1, NUM_OF_SPHERES, canonized);

	double newPoints[3*NUM_OF_SPHERES];
	for(int i=0; i< NUM_OF_SPHERES; i++){
		for(int j=0; j<3; j++){
			newPoints[3*i+j] = (this->p) (3*lab[i]+j);
		}
	}
	memcpy(&(this->p), newPoints, sizeof(double)*NUM_OF_SPHERES*3);
	memcpy(this->g, canonized, sizeof(graph)*NUM_OF_SPHERES);
	
	densenauty(this->g, lab, ptn, orbits, &options, &stats, 1, NUM_OF_SPHERES, canonized);
	
	group = groupptr(FALSE);
	makecosetreps(group);
	ptype.np=0;
	allgroup3(group,countperms,&ptype);

	
	
	
	chooseTriangle();
	return fixTriangle();
	


	
}
