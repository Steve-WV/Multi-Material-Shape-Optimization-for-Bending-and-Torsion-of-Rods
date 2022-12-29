#include <unistd.h>
#include <iostream>
#include "mesh.hpp"
#include "vertex.hpp"
#include "tri.hpp"
#include "pde.hpp"
#include "tmv.hpp"

#include <Eigen/Dense>



int main(int argc, char *argv[]) {
	
	std::string obj_fname, opt_fname="";
	bool vis=false;
	int c;
	while ((c = getopt (argc, argv, "f:c:v")) != -1)
		switch (c)
	{
		case 'f':
		obj_fname = std::string(optarg);
		case 'c':
		opt_fname = std::string(optarg);
		break;
		case 'v':
		vis = true;
		break;
		case '?':
		if (optopt == 'f' || optopt == 'c')
			std::cerr << "Option -" << optopt <<" requires an argument.\n";
		else if (isprint (optopt))
			std::cerr << "Unknown option `-" << optopt << "'.\n";
		else
			std::cerr << "Unknown option character `" << optopt << "'.\n";
		return 1;
		default:
		return(-1);
	}
    
	if (argc-optind != 0 ||  obj_fname == "") {
		std::cerr << "Problem with options.";
		return(-1);
	}
		
		
	Mesh* m = new Mesh();
	m->Load(obj_fname);
	ImplicitPDESystem* pde = new ImplicitPDESystem(m);
	pde->Init(opt_fname, vis);
	
	pde->Solve(); 
	
}
