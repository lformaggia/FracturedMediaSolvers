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

		gmm::size_type flag(0);
		getline (myfile,line);
		getline (myfile,line);
		getline (myfile,line);

		M_isMetric=false;
		gmm::size_type pos11;
		pos11=line.find_first_of("=");
	        char provv11[1];
		gmm::size_type length11=line.copy(provv11, 1, pos11+2);

		if (provv11=="M") {M_isMetric=true; std::cout << "METRICO"<<std::endl;}
		while(myfile.good() && flag==0){
	
			getline (myfile,line);
			gmm::size_type pos, pos1;
			pos=line.find("No_Fractures");
 			if (pos!=std::string::npos){
			pos=line.find_first_of("0123456789");
			pos1=line.find_last_of("0123456789");
			char provv[pos1-pos];
			gmm::size_type length=line.copy(provv, pos1+1-pos , pos);
			provv[length]='\0';
			M_nfractures=atoi(provv);
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
				char provv[found-found2];
				gmm::size_type length=line.copy(provv, found-found2 , found2+1);
				provv[length]='\0';
				 pz=atof(provv);
				line.resize(found2+1);

				found=line.find_last_of("0123456789");
				line.resize(found+1);
				found2=line.find_last_of(" ");
				char provv1[found-found2];
				length=line.copy(provv1, found-found2 , found2+1);
				provv1[length]='\0';
				py=atof(provv1);
				line.resize(found2+1);

				found=line.find_last_of("0123456789");
				line.resize(found+1);
				found2=line.find_last_of(" ");
				char provv2[found-found2];
				length=line.copy(provv2, found-found2 , found2+1);
				provv2[length]='\0';
				px=atof(provv2);
}

void readProperties(std::string line, gmm::size_type &n_points, Real &perm, Real  &compr, Real &aperture){
	gmm::size_type found=line.find_last_of("0123456789");
				line.resize(found+1);
				gmm::size_type found2=line.find_last_of(" ");
				char provv[found-found2];
				gmm::size_type length=line.copy(provv, found-found2 , found2+1);
				provv[length]='\0';
				aperture=atof(provv);
				line.resize(found2+1);

				found=line.find_last_of("0123456789");
				line.resize(found+1);
				found2=line.find_last_of(" ");
				char provv1[found-found2];
				length=line.copy(provv1, found-found2 , found2+1);
				provv1[length]='\0';
				compr=atof(provv1);
				line.resize(found2+1);

				found=line.find_last_of("0123456789");
				line.resize(found+1);
				found2=line.find_last_of(" ");
				char provv2[found-found2];
				length=line.copy(provv2, found-found2 , found2+1);
				provv2[length]='\0';
				perm=atof(provv2);
				line.resize(found2+1);

				found=line.find_last_of("0123456789");
				line.resize(found+1);
				found2=line.find_last_of(" ");
				line.resize(found2+1);
				
				found=line.find_last_of("0123456789");
				line.resize(found+1);
				found2=line.find_last_of(" ");
				char provv3[found-found2];
				length=line.copy(provv3, found-found2 , found2+1);
				provv3[length]='\0';
				n_points=atoi(provv3);
	
				line.resize(found2+1);
}

} // namespace Geometry
