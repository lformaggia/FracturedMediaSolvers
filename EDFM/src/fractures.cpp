 /*!
 *	@file geomFault.cpp
 *	@brief Class for Fault (definition).
 *
 *	@author Luca Turconi <lturconi@gmail.com>
 *  @date 22-09-2012
 *
 */

#include<omp.h>
//#include "/usr/lib/gcc/x86_64-unknown-linux-gnu/4.7.1/include/omp.h"

#include "fractures.hpp"


namespace Geometry
{
	Fractures::Fractures(const std::string nomefile){
		std::cout << "sono qui"<<std::endl;
		std::ifstream myfile;
		std::string line;
		std::string nomefile1("data/");
		nomefile1.append(nomefile);

		myfile.open(nomefile1.c_str());
		if(myfile.fail())
		  {
		    std::cerr<<"Error opening file "<<nomefile1<<std::endl;
		    std::exit(1);
		  }

		gmm::size_type flag(0);
		getline (myfile,line);
		getline (myfile,line);
		getline (myfile,line);

		M_isMetric=false;
		gmm::size_type pos11;
		pos11=line.find_first_of("=");
		std::string provv11=line.substr(pos11+2,1);
	    //    char provv11[1];
		//line.copy(provv11, 1, pos11+2);

		if (provv11==std::string("M")) {M_isMetric=true; std::cout << "METRICO"<<std::endl;}
		while(myfile.good() && flag==0){

			getline (myfile,line);
			gmm::size_type pos, pos1;
			pos=line.find("No_Fractures");
 			if (pos!=std::string::npos){
			pos=line.find_first_of("0123456789");
			pos1=line.find_last_of("0123456789");
			std::string provv=line.substr(pos,pos1+1-pos);
			M_nfractures=atoi(provv.c_str());
			}
			pos=line.find("BEGIN FRACTURE");
 			if (pos!=std::string::npos){

			flag=1;
			}
		}

		for (gmm::size_type i=0; i<M_nfractures;++i ){

			getline (myfile,line);
			gmm::size_type n_punti;
			Real perm, compr, aperture;
			readProperties(line, n_punti, perm, compr, aperture);
			std::vector<Point3D> punti;
			for (gmm::size_type j=0; j<n_punti;++j){
				getline (myfile,line);
				Real px,py,pz;
				readPoints(line,px,py,pz);
				Point3D punto(px,py,pz);
				punti.push_back(punto);
				//std::cout << punto<<std::endl;
			}
			getline (myfile,line);

			if (n_punti==3) {
			Point3D PP(0.55*punti[punti.size()-1]+0.55*punti[punti.size()-2]-0.1*punti[punti.size()-3]);
			Point3D D(punti[punti.size()-1]);
			punti.pop_back();
			punti.push_back(PP);
			punti.push_back(D);
			}
			if (n_punti==5) {punti.pop_back();}

			Fracture frattura(punti[0], punti[1], punti[2], punti[3]);
			frattura.setMetric(M_isMetric);
			M_fractures.push_back(frattura);
			M_fractures[i].setProperties(perm, compr, aperture);

		}

		myfile.clear();              // forget we hit the end of file
		myfile.seekg(0, std::ios::beg);   // move to the start of the file
	}

	void Fractures::computeIntersections(bool completo)
	{

	for (gmm::size_type i=0; i<M_nfractures;++i){

		for (gmm::size_type j=0; j<M_nfractures;++j && j!=i){
		std::vector<Point3D> puntiint;
		bool primo(M_fractures[i].isIntersectedBy(M_fractures[j], puntiint));
		bool secondo( M_fractures[j].isIntersectedBy(M_fractures[i], puntiint));

			if (primo||secondo )
					{
						M_fractures[i].setInt(j, M_fractures[j], puntiint, completo);
					}
		}
	}
	}




//------------------distruttore-------------------------

	Fractures::~Fractures() {}


void readPoints(std::string line, Real &px, Real &py, Real &pz){
	gmm::size_type found=line.find_last_of("0123456789");
				line.resize(found+1);
				gmm::size_type found2=line.find_last_of(" ");
				std::string provv=line.substr(found2+1,found-found2);
				pz=atof(provv.c_str());
				line.resize(found2+1);

				found=line.find_last_of("0123456789");
				line.resize(found+1);
				found2=line.find_last_of(" ");
				std::string provv1=line.substr(found2+1,found-found2);
				py=atof(provv1.c_str());
				line.resize(found2+1);

				found=line.find_last_of("0123456789");
				line.resize(found+1);
				found2=line.find_last_of(" ");
				std::string provv2=line.substr(found2+1,found-found2);
				px=atof(provv2.c_str());
}

void readProperties(std::string line, gmm::size_type &n_points, Real &perm, Real  &compr, Real &aperture){
	gmm::size_type found=line.find_last_of("0123456789");
				line.resize(found+1);
				gmm::size_type found2=line.find_last_of(" ");
				std::string provv=line.substr(found2+1,found-found2);
				aperture=atof(provv.c_str());
				line.resize(found2+1);

				found=line.find_last_of("0123456789");
				line.resize(found+1);
				found2=line.find_last_of(" ");
				std::string provv1=line.substr(found2+1,found-found2);
				compr=atof(provv1.c_str());
				line.resize(found2+1);

				found=line.find_last_of("0123456789");
				line.resize(found+1);
				found2=line.find_last_of(" ");
				std::string provv2=line.substr(found2+1,found-found2);
				perm=atof(provv2.c_str());
				line.resize(found2+1);

				found=line.find_last_of("0123456789");
				line.resize(found+1);
				found2=line.find_last_of(" ");
				line.resize(found2+1);

				found=line.find_last_of("0123456789");
				line.resize(found+1);
				found2=line.find_last_of(" ");
				provv2=line.substr(found2+1,found-found2);
				n_points=atoi(provv2.c_str());

				line.resize(found2+1);
}

} // namespace Geometry
