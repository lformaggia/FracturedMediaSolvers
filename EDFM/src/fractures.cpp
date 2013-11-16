/*!
*  @file geomFault.cpp
*  @brief Class for Fault (definition).
*
*  @author Luca Turconi <lturconi@gmail.com>
*  @date 22-09-2012
*
*/

#include<omp.h>
//#include "/usr/lib/gcc/x86_64-unknown-linux-gnu/4.7.1/include/omp.h"

#include "fractures.hpp"


namespace Geometry
{
  Fractures::Fractures (const std::string nomefile, Real conv_z, std::string direzione, Real angle)
  {

    std::ifstream myfile;
    std::string line;
    //std::string nomefile1("data/");
    //nomefile1.append(nomefile);

    myfile.open (nomefile.c_str() );
    if (myfile.fail() )
    {
      std::cerr << "Error opening file " << nomefile << std::endl;
      std::exit (1);
    }

    gmm::size_type flag (0);
    getline (myfile, line);
    getline (myfile, line);
    getline (myfile, line);

    M_isMetric = false;
    gmm::size_type pos11;
    pos11 = line.find_first_of ("=");
    char provv11[1];
    line.copy (provv11, 1, pos11 + 2);

    if (provv11[0] == 'M')
    {
      M_isMetric = true;
      std::cout << "METRICO" << std::endl;
    }
    while (myfile.good() && flag == 0)
    {

      getline (myfile, line);
      gmm::size_type pos, pos1;
      pos = line.find ("No_Fractures");
      if (pos != std::string::npos)
      {
        pos = line.find_first_of ("0123456789");
        pos1 = line.find_last_of ("0123456789");
	//        char provv[pos1 - pos];
        //gmm::size_type length = line.copy (provv, pos1 + 1 - pos , pos);
        //provv[length] = '\0';
	std::string temp=line.substr(pos,pos1+1-pos);
        //M_nfractures = atoi (provv);
        M_nfractures = atoi (temp.c_str());
      }
      pos = line.find ("BEGIN FRACTURE");
      if (pos != std::string::npos)
      {

        flag = 1;
      }
    }

    for (gmm::size_type i = 0; i < M_nfractures; ++i )
    {

      getline (myfile, line);
      gmm::size_type n_punti;
      Real perm, compr, aperture;
      readProperties (line, n_punti, perm, compr, aperture);
      std::vector<Point3D> punti;
      for (gmm::size_type j = 0; j < n_punti; ++j)
      {
        getline (myfile, line);
        Real px, py, pz;
        readPoints (line, px, py, pz);
        Point3D punto (px, py, conv_z * pz);
        punti.push_back (punto);

      }
      getline (myfile, line);

      if (n_punti == 3)
      {
        Point3D PP (0.55 * punti[punti.size() - 1] + 0.55 * punti[punti.size() - 2] - 0.1 * punti[punti.size() - 3]);
        Point3D D (punti[punti.size() - 1]);
        punti.pop_back();
        punti.push_back (PP);
        punti.push_back (D);
      }
      if (n_punti == 5)
      {
        punti.pop_back();
      }
      Point3D centro (0.25 * punti[0] + 0.25 * punti[1] + 0.25 * punti[2] + 0.25 * punti[3]);

      Fracture frattura (apply_shear (punti[0] - centro, angle, direzione) + centro, apply_shear (punti[1] - centro, angle, direzione) + centro, apply_shear (punti[2] - centro, angle, direzione) + centro, apply_shear (punti[3] - centro, angle, direzione) + centro);
      //punti[1]=punti[1]+punti[0];
      //punti[2]=punti[2]+punti[0];
      //punti[3]=punti[3]+punti[0];
      frattura.setMetric (M_isMetric);
      M_fractures.push_back (frattura);
      M_fractures[i].setProperties (perm, compr, aperture);

    }

    myfile.clear();              // forget we hit the end of file
    myfile.seekg (0, std::ios::beg);  // move to the start of the file
  }

  void Fractures::computeIntersections (bool completo)
  {

    for (gmm::size_type i = 0; i < M_nfractures; ++i)
    {

      for (gmm::size_type j = 0; j < M_nfractures; ++j)
      {
        std::vector<Point3D> puntiint;
        bool primo (M_fractures[i].isIntersectedBy (M_fractures[j], puntiint) );
        bool secondo ( M_fractures[j].isIntersectedBy (M_fractures[i], puntiint) );
        if (puntiint.size() > 0)
        {
          Segment SS (maxSegment (puntiint) );

          if ( (primo || secondo) && (i != j) && (M_fractures[i].areaFault() > 0) && (M_fractures[j].areaFault() > 0) && SS.length() > 	                		 		       EDFM_Tolerances::POINT2DCOINCIDENT)
          {
            M_fractures[i].setInt (j, M_fractures[j], puntiint, completo);

          }
        }
      }

    }
    if (completo == true)
    {
      this->setInterFTransm(); //calcolo delle medie armoniche.
    }
  }




  //------------------distruttore-------------------------

  Fractures::~Fractures() {}


  void readPoints (std::string line, Real& px, Real& py, Real& pz)
  {
    gmm::size_type found = line.find_last_of ("0123456789");
    line.resize (found + 1);
    gmm::size_type found2 = line.find_last_of (" ");
    std::string temp=line.substr(found2+1,found-found2);
    //char provv[found - found2];
    //    gmm::size_type length = line.copy (provv, found - found2 , found2 + 1);
    //provv[length] = '\0';
    //pz = atof (provv);
    pz = atof (temp.c_str());
    line.resize (found2 + 1);

    found = line.find_last_of ("0123456789");
    line.resize (found + 1);
    found2 = line.find_last_of (" ");
    //char provv1[found - found2];
    //length = line.copy (provv1, found - found2 , found2 + 1);
    temp=line.substr(found2+1,found-found2);
    //    provv1[length] = '\0';
    //py = atof (provv1);
    py = atof (temp.c_str());
    line.resize (found2 + 1);

    found = line.find_last_of ("0123456789");
    line.resize (found + 1);
    found2 = line.find_last_of (" ");
    temp=line.substr(found2+1,found-found2);
    //char provv2[found - found2];
    // length = line.copy (provv2, found - found2 , found2 + 1);
    //provv2[length] = '\0';
    //px = atof (provv2);
    px=atof(temp.c_str());
  }

  void readProperties (std::string line, gmm::size_type& n_points, Real& perm, Real&  compr, Real& aperture)
  {
    gmm::size_type found = line.find_last_of ("0123456789");
    line.resize (found + 1);
    gmm::size_type found2 = line.find_last_of (" ");
    std::string temp=line.substr(found2+1,found-found2);
    //char provv[found - found2];
    //gmm::size_type length = line.copy (provv, found - found2 , found2 + 1);
    //provv[length] = '\0';
    //aperture = atof (provv);
    aperture=atof(temp.c_str());
    line.resize (found2 + 1);

    found = line.find_last_of ("0123456789");
    line.resize (found + 1);
    found2 = line.find_last_of (" ");
    temp=line.substr(found2+1,found-found2);
    //char provv1[found - found2];
    //length = line.copy (provv1, found - found2 , found2 + 1);
    //provv1[length] = '\0';
    //compr = atof (provv1);
    compr = atof (temp.c_str());
    line.resize (found2 + 1);

    found = line.find_last_of ("0123456789");
    line.resize (found + 1);
    found2 = line.find_last_of (" ");
    temp  = line.substr(found2+1,found-found2);
    //char provv2[found - found2];
    //length = line.copy (provv2, found - found2 , found2 + 1);
    //provv2[length] = '\0';
    //perm = atof (provv2);
    perm = atof (temp.c_str());
    
    line.resize (found2 + 1);

    found = line.find_last_of ("0123456789");
    line.resize (found + 1);
    found2 = line.find_last_of (" ");
    line.resize (found2 + 1);

    found = line.find_last_of ("0123456789");
    line.resize (found + 1);
    found2 = line.find_last_of (" ");
    temp =  line.substr(found2+1,found-found2);
    //char provv3[found - found2];
    //length = line.copy (provv3, found - found2 , found2 + 1);
    //provv3[length] = '\0';
    //n_points = atoi (provv3);
    n_points = atoi (temp.c_str());

    line.resize (found2 + 1);
  }

  void Fractures::setInterFTransm()
  {
    //ciclo sulle fratture
    for (gmm::size_type i = 0; i < M_nfractures; ++i)
    {
      std::vector<gmm::size_type> quali (M_fractures[i].getIsInt() );
      for (gmm::size_type j = 0; j < quali.size(); ++j)
      {
        std::vector<gmm::size_type> quali_altra (M_fractures[quali[j]].getIsInt() );
        std::vector<gmm::size_type>::iterator it;
        it = find (quali_altra.begin(), quali_altra.end(), i);
        gmm::size_type jj (it - quali_altra.begin() );

        for (gmm::size_type k = 0; k < M_fractures[i].inter() [j].M_i.size(); ++k)
        {
          for (gmm::size_type h = 0; h < M_fractures[quali[j]].inter() [jj].M_i.size(); ++h)
          {
            if (M_fractures[i].inter() [j].Which2[k] == M_fractures[quali[j]].inter() [jj].Which1[h])
            {
              Real T1, T2, Tave;
              T1 = M_fractures[i].inter() [j].M_Tf1f2h[k];
              T2 = M_fractures[quali[j]].inter() [jj].M_Tf1f2h[h];
              Tave = 1. / (1. / T1 + 1. / T2);
              M_fractures[i].setInterTransm (j, k, Tave);
              M_fractures[quali[j]].setInterTransm (jj, h, Tave);
            }
          }
        }




      }//ciclo su quelle che interseca
    }//fine ciclo sulle fratture

  }

  void Fractures::exportIntersectionMatrix (std::ofstream& myfile)
  {
    for (gmm::size_type i = 0; i < M_nfractures; ++i)
    {
      for (gmm::size_type j = 0; j < M_nfractures; ++j)
      {
        std::vector<gmm::size_type> intersezioni (this->M_fractures[i].getIsInt() );
        if (find (intersezioni.begin(), intersezioni.end(), j) != intersezioni.end() )
        {
          myfile << 1 << "\t";
        }
        else
        {
          myfile << 0 << "\t";
        }
      }
      myfile << "\n";
    }
  }

} // namespace Geometry
