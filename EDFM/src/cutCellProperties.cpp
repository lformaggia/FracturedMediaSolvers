/*!
*  @file geomTetra.cpp
*  @brief Class for Tetra in 3D space (definition).
*
*
*/

#include<cmath>
#include<limits>
#include<fstream>
#include<iomanip>
#include<vector>
#include<algorithm>
#include "cutCellProperties.hpp"


namespace Geometry
{

  // --------------------   Class CProp --------------------

  // ==================================================
  // Constructors & Destructor
  // ==================================================
  CProp::CProp() : M_vol(), M_aree(), M_dmedio(), M_iteratorcellsbegin(), M_iteratorcellsend(), M_gridpointer(), M_faultpointer() {}

  //------------------costruttore


  CProp::CProp (const Intersect::GridIntersections& storage, CPgrid* gridpointer, Fracture* faultpointer)
  {
    M_iteratorcellsbegin = storage.begin();
    M_iteratorcellsend = storage.end();
    M_gridpointer = gridpointer;
    M_faultpointer = faultpointer;

    M_vol.resize (M_gridpointer->Nx() *M_gridpointer->Ny() *M_gridpointer->Nz() );
    M_aree.resize (M_gridpointer->Nx() *M_gridpointer->Ny() *M_gridpointer->Nz() );
    M_CG.resize (M_gridpointer->Nx() *M_gridpointer->Ny() *M_gridpointer->Nz() );
    M_dmedio.resize (M_gridpointer->Nx() *M_gridpointer->Ny() *M_gridpointer->Nz() );

  }


  //------------------distruttore-------------------------

  CProp::~CProp() {}

  //--------------settare le proprietà--------------------

  void CProp::setProperties()
  {

    M_Ne = 0;
    gmm::size_type counter (0);
    for (Intersect::GridIntersections_Const_Iterator_Type it = M_iteratorcellsbegin;
         it != M_iteratorcellsend; ++it)
    {
      counter += 1;
    }

    for (gmm::size_type i = 0; i < M_faultpointer->getIsInt().size(); ++i)
    {

      std::vector<Real> ooo (counter, 0.);
      M_dmedioint.push_back (ooo);
    }

    // Loop on cells interested by the intersection
    for (Intersect::GridIntersections_Const_Iterator_Type it = M_iteratorcellsbegin;
         it != M_iteratorcellsend; ++it)
    {
      M_Ne = M_Ne + 1;

      M_i.push_back (it->second.i() );
      M_j.push_back (it->second.j() );
      M_k.push_back (it->second.k() );
      // Create the cell
      CPcell cella (M_gridpointer->cell ( (*it).second.i(), (*it).second.j(), (*it).second.k() ) );
      // Intersection points: true and virtual
      std::vector<Point3D> punti1 (this->getIntPoints (it) );
      // True vs virtial intersection points.
      std::vector<bool> puntiIsReal (this->getIsIntPointReal (it) );
      
      gmm::size_type npunti_faglia (punti1.size() );

      Point3D Aa (punti1[0]);
      // Compute baricenter of all intersection points
      for (gmm::size_type i = 1; i < punti1.size(); ++i)
      {
        Aa = Aa + punti1[i];
      }
      Aa.x = Aa.x / Real (npunti_faglia);
      Aa.y = Aa.y / Real (npunti_faglia);
      Aa.z = Aa.z / Real (npunti_faglia);

      // Compute normal
      Point3D no (M_faultpointer->normal (0.5, 0.5) );
      Real typicalL;
      // A reference length
      typicalL = (cella.getVertex (1) - cella.getVertex (2) ).norm();

      std::vector<Real> uv (M_faultpointer->inv_param (Aa) );

      std::vector<Point3D> punti2 (punti1);

      std::vector<Point3D> puntiarea (punti1);
      std::vector<Point3D> puntiAreaNew;
      std::vector<bool>::iterator location;
      std::vector<bool>::iterator  location2;

      // Here we partition "real" intersection points with "virtual" ones
      // Of course the idea is that thet are consecutively grouped into
      // puntiIsReal
      location = std::find ( puntiIsReal.begin(), puntiIsReal.end(), false);
      location2 = std::find ( puntiIsReal.begin(), puntiIsReal.end(), true);

      buildIntSegments (it);

      bool isPartial (false);
      if (true) //(location2!=puntiIsReal.end()){
      {

        if (true) //(location!=puntiIsReal.end()){
        {
          puntiAreaNew = this->addPoints4area (puntiarea, puntiIsReal, it); //il pizzino
          isPartial = true;

        }
        else
        {
          puntiAreaNew = puntiarea;
        }
        // computing baricenter as average.
        M_CG[ (*it).first] = this->setCG (puntiAreaNew);

        M_puntiAree.push_back (puntiAreaNew);

	// Add a off-plane point to generate 3D convex hull with QHull
        puntiAreaNew.push_back (this->setCG (puntiAreaNew) + no * typicalL);
        // WRAPPER QHULL
        Hull faccia (puntiAreaNew);
	// ! For debugging
	//        if (it->second.i() == 62 && it->second.j() == 27 && it->second.k() == 3)
        //{
	// std::cout << npunti_faglia << std::endl;
        //}
        if (faccia.getNtetra() > 0)
        {
	  // Computes area
          M_aree[ (*it).first] = this->setIntArea (faccia, puntiAreaNew.size() - 1);

          if (location != puntiIsReal.end() )
          {
	    // Computes barycenter
            M_CG[ (*it).first] = this->setCG (faccia, puntiAreaNew.size() - 1);
          }
        }

	// Take out last point (could have been done earlier 
        puntiAreaNew.pop_back();
        for (gmm::size_type i = 0; i < M_faultpointer->getIsInt().size(); ++i)
        {

          Real dmediosingolo (0);

          if (faccia.getNtetra() > 0 )

          {
	    // This beats me! 
            dmediosingolo = this->setIntdist_linea (puntiAreaNew, no, M_faultpointer->inter() [i], it, isPartial);

          }

          M_dmedioint[i][M_Ne - 1] = dmediosingolo;

        }

      }
      Aa = M_CG[ (*it).first];
      this->getCellPoints (it, punti1, 1, Aa);

      Hull calimero (punti1);
      if (calimero.getNtetra() > 0)
      {
        M_vol[ (*it).first] += calimero.getVolume();
        M_dmedio[ (*it).first] += this->setIntd (calimero, Aa);
      }

      this->getCellPoints (it, punti2, -1, Aa);

      Hull calimero2 (punti2);

      if (calimero2.getNtetra() > 0)
      {
        M_vol[ (*it).first] += calimero2.getVolume();
        M_dmedio[ (*it).first] += this->setIntd (calimero2, Aa);
      }

      if (M_vol[ (*it).first] > 0)
      {
        M_dmedio[ (*it).first] = M_dmedio[ (*it).first] / M_vol[ (*it).first];
      }
      else
      {
        M_dmedio[ (*it).first] = 0;
      }

    }



  }

  //--------------segmenti di intersezione........................................

  void CProp::buildIntSegments (Intersect::GridIntersections_Const_Iterator_Type& it)
  {
    std::vector<Point3D> puntiX1, puntiY1, puntiZ1;
    std::vector<bool> isXreal1, isYreal1, isZreal1;
    std::vector<Point3D> puntiX2, puntiY2, puntiZ2;
    std::vector<bool> isXreal2, isYreal2, isZreal2;
    std::vector<Segment> latiFaglia;
    Segment ss1 (M_faultpointer->A(), M_faultpointer->B() );
    latiFaglia.push_back (ss1);
    Segment ss2 (M_faultpointer->B(), M_faultpointer->C() );
    latiFaglia.push_back (ss2);
    Segment ss3 (M_faultpointer->C(), M_faultpointer->D() );
    latiFaglia.push_back (ss3);
    Segment ss4 (M_faultpointer->D(), M_faultpointer->A() );
    latiFaglia.push_back (ss4);
    CPcell cella (M_gridpointer->cell ( (*it).second.i(), (*it).second.j(), (*it).second.k() ) );

    Point3D PP (0, 0, 0);
    Segment fake (PP, PP);

    //--------------------faccia X1----------------------------------------

    for (Intersect::CellIntersections_Const_Iterator_Type jt = (*it).second.begin();
         jt != (*it).second.end(); ++jt)
    {
      if (jt->first == 4 || jt->first == 8 || jt->first == 5 || jt->first == 12)
      {
        puntiX1.push_back (jt->second);
        isXreal1.push_back (true);
      }
      if (jt->first == 104 || jt->first == 108 || jt->first == 105 || jt->first == 112)
      {
        puntiX1.push_back (jt->second);
        isXreal1.push_back (false);
      }
    }
    if (puntiX1.size() == 2)
    {
      Segment provv (puntiX1[0], puntiX1[1]);
      for (gmm::size_type i = 0; i < 4; ++i)
      {
        bool giafatto (false);
        Point3D comodo;
        if (false) // (isXreal1[0]==false && isXreal1[1]==false)
        {
          puntiX1[0] = PP;
          puntiX1[1] = PP;
        }
        else
        {
          if (isXreal1[0] == false)
          {
            if (latiFaglia[i].intersectTheSegment (provv, comodo) && provv.isIn (comodo) )
            {
              giafatto = true;
              puntiX1[0] = comodo;
            }
          }
          if (isXreal1[1] == false && giafatto == false)
          {
            if (latiFaglia[i].intersectTheSegment (provv, comodo) && provv.isIn (comodo) )
            {
              puntiX1[1] = comodo;

            }
          }
        }
      }
      Segment ss (puntiX1[0], puntiX1[1]);
      M_S1x.push_back (ss);
    }
    else
    {
      M_S1x.push_back (fake);
    }

    //--------------------faccia Y1----------------------------------------

    for (Intersect::CellIntersections_Const_Iterator_Type jt = (*it).second.begin();
         jt != (*it).second.end(); ++jt)
    {
      if (jt->first == 1 || jt->first == 5 || jt->first == 6 || jt->first == 9)
      {
        puntiY1.push_back (jt->second);
        isYreal1.push_back (true);
      }
      if (jt->first == 101 || jt->first == 105 || jt->first == 106 || jt->first == 109)
      {
        puntiY1.push_back (jt->second);
        isYreal1.push_back (false);
      }
    }
    if (puntiY1.size() == 2)
    {
      Segment provv (puntiY1[0], puntiY1[1]);
      for (gmm::size_type i = 0; i < 4; ++i)
      {
        bool giafatto (false);
        Point3D comodo;
        if (false) //(isYreal1[0]==false && isYreal1[1]==false)
        {
          puntiY1[0] = PP;
          puntiY1[1] = PP;
        }
        else
        {
          if (isYreal1[0] == false)
          {
            if (latiFaglia[i].intersectTheSegment (provv, comodo) && provv.isIn (comodo) )
            {
              giafatto = true;
              puntiY1[0] = comodo;
            }
          }
          if (isYreal1[1] == false && giafatto == false)
          {
            if (latiFaglia[i].intersectTheSegment (provv, comodo) && provv.isIn (comodo) )
            {
              puntiY1[1] = comodo;
            }
          }
        }
      }
      Segment ss (puntiY1[0], puntiY1[1]);
      M_S1y.push_back (ss);
    }
    else
    {
      M_S1y.push_back (fake);
    }

    //--------------------faccia Z1----------------------------------------

    for (Intersect::CellIntersections_Const_Iterator_Type jt = (*it).second.begin();
         jt != (*it).second.end(); ++jt)
    {

      if (jt->first == 1 || jt->first == 2 || jt->first == 3 || jt->first == 4)
      {
        puntiZ1.push_back (jt->second);
        isZreal1.push_back (true);
      }
      if (jt->first == 101 || jt->first == 102 || jt->first == 103 || jt->first == 104)
      {
        puntiZ1.push_back (jt->second);
        isZreal1.push_back (false);
      }
    }

    if (puntiZ1.size() == 2)
    {
      Segment provv (puntiZ1[0], puntiZ1[1]);
      gmm::size_type cont (0);
      for (gmm::size_type i = 0; i < 4; ++i)
      {
        bool giafatto (false);
        Point3D comodo;
        if (isZreal1[0] == false && isZreal1[1] == false && cont < 2)
        {
          //puntiZ1[0]=PP;
          //puntiZ1[1]=PP;
          if (latiFaglia[i].intersectTheSegment (provv, comodo) && provv.isIn (comodo) )
          {
            puntiZ1[cont] = comodo;
            cont = cont + 1;
          }
        }
        else
        {
          if (isZreal1[0] == false)
          {
            if (latiFaglia[i].intersectTheSegment (provv, comodo) && provv.isIn (comodo) )
            {
              giafatto = true;
              puntiZ1[0] = comodo;
            }
          }
          if (isZreal1[1] == false && giafatto == false)
          {
            if (latiFaglia[i].intersectTheSegment (provv, comodo) && provv.isIn (comodo) )
            {
              puntiZ1[1] = comodo;
            }
          }
        }
      }
      Segment ss (puntiZ1[0], puntiZ1[1]);
      M_S1z.push_back (ss);
    }
    else
    {
      M_S1z.push_back (fake);
    }

    //-------------------------faccia X2-----------------------------------------

    for (Intersect::CellIntersections_Const_Iterator_Type jt = (*it).second.begin();
         jt != (*it).second.end(); ++jt)
    {
      if (jt->first == 6 || jt->first == 7 || jt->first == 10 || jt->first == 2)
      {
        puntiX2.push_back (jt->second);
        isXreal2.push_back (true);
      }
      if (jt->first == 106 || jt->first == 107 || jt->first == 110 || jt->first == 102)
      {
        puntiX2.push_back (jt->second);
        isXreal2.push_back (false);
      }
    }
    if (puntiX2.size() == 2)
    {
      Segment provv (puntiX2[0], puntiX2[1]);
      for (gmm::size_type i = 0; i < 4; ++i)
      {
        bool giafatto (false);
        Point3D comodo;
        if (false) //(isXreal2[0]==false && isXreal2[1]==false)
        {
          puntiX2[0] = PP;
          puntiX2[1] = PP;
        }
        else
        {
          if (isXreal2[0] == false)
          {
            if (latiFaglia[i].intersectTheSegment (provv, comodo) && provv.isIn (comodo) )
            {
              giafatto = true;
              puntiX2[0] = comodo;
            }
          }
          if (isXreal2[1] == false && giafatto == false)
          {
            if (latiFaglia[i].intersectTheSegment (provv, comodo) && provv.isIn (comodo) )
            {
              puntiX2[1] = comodo;

            }
          }
        }
      }
      Segment ss (puntiX2[0], puntiX2[1]);
      M_S2x.push_back (ss);
    }
    else
    {
      M_S2x.push_back (fake);
    }
    //---------------------------faccia Y2----------------------------------------------

    for (Intersect::CellIntersections_Const_Iterator_Type jt = (*it).second.begin();
         jt != (*it).second.end(); ++jt)
    {
      if (jt->first == 3 || jt->first == 7 || jt->first == 8 || jt->first == 11)
      {
        puntiY2.push_back (jt->second);
        isYreal2.push_back (true);
      }
      if (jt->first == 103 || jt->first == 107 || jt->first == 108 || jt->first == 111)
      {
        puntiY2.push_back (jt->second);
        isYreal2.push_back (false);
      }
    }
    if (puntiY2.size() == 2)
    {
      Segment provv (puntiY2[0], puntiY2[1]);
      for (gmm::size_type i = 0; i < 4; ++i)
      {
        bool giafatto (false);
        Point3D comodo;
        if (false) //(isYreal2[0]==false && isYreal2[1]==false)
        {
          puntiY2[0] = PP;
          puntiY2[1] = PP;
        }
        else
        {
          if (isYreal2[0] == false)
          {
            if (latiFaglia[i].intersectTheSegment (provv, comodo) && provv.isIn (comodo) )
            {
              giafatto = true;
              puntiY2[0] = comodo;
            }
          }
          if (isYreal2[1] == false && giafatto == false)
          {
            if (latiFaglia[i].intersectTheSegment (provv, comodo) && provv.isIn (comodo) )
            {
              puntiY2[1] = comodo;
            }
          }
        }
      }
      Segment ss (puntiY2[0], puntiY2[1]);
      M_S2y.push_back (ss);
    }
    else
    {
      M_S2y.push_back (fake);
    }

    //--------------------faccia Z2----------------------------------------

    for (Intersect::CellIntersections_Const_Iterator_Type jt = (*it).second.begin();
         jt != (*it).second.end(); ++jt)
    {

      if (jt->first == 9 || jt->first == 10 || jt->first == 11 || jt->first == 12)
      {
        puntiZ2.push_back (jt->second);
        isZreal2.push_back (true);
      }
      if (jt->first == 109 || jt->first == 110 || jt->first == 111 || jt->first == 112)
      {
        puntiZ2.push_back (jt->second);
        isZreal2.push_back (false);
      }
    }

    if (puntiZ2.size() == 2)
    {
      Segment provv (puntiZ2[0], puntiZ2[1]);
      gmm::size_type cont (0);
      for (gmm::size_type i = 0; i < 4; ++i)
      {
        bool giafatto (false);
        Point3D comodo;
        if  (isZreal2[0] == false && isZreal2[1] == false && cont < 2)
        {
          //puntiZ2[0]=PP;
          //puntiZ2[1]=PP;
          if (latiFaglia[i].intersectTheSegment (provv, comodo) && provv.isIn (comodo) )
          {
            puntiZ2[cont] = comodo;
            cont = cont + 1;
          }
        }
        else
        {
          if (isZreal2[0] == false)
          {
            if (latiFaglia[i].intersectTheSegment (provv, comodo) && provv.isIn (comodo) )
            {
              giafatto = true;
              puntiZ2[0] = comodo;
            }
          }
          if (isZreal2[1] == false && giafatto == false)
          {
            if (latiFaglia[i].intersectTheSegment (provv, comodo) && provv.isIn (comodo) )
            {
              puntiZ2[1] = comodo;
            }
          }
        }
      }
      Segment ss (puntiZ2[0], puntiZ2[1]);
      M_S2z.push_back (ss);
    }
    else
    {
      M_S2z.push_back (fake);
    }

  }

  //--------------prendere per ogni cella i punti di intersezione-------------------

  std::vector<Point3D> CProp::getIntPoints (Intersect::GridIntersections_Const_Iterator_Type& it)
  {

    std::vector<Point3D> puntiInt;
    for (Intersect::CellIntersections_Const_Iterator_Type jt = (*it).second.begin();
         jt != (*it).second.end(); ++jt)
    {
      puntiInt.push_back (jt->second);
    }

    return puntiInt;
  }

  //--------------prendere per ogni cella i punti di intersezione-------------------

  std::vector<bool> CProp::getIsIntPointReal (Intersect::GridIntersections_Const_Iterator_Type& it)
  {

    std::vector<bool> IsReal;
    for (Intersect::CellIntersections_Const_Iterator_Type jt = (*it).second.begin();
         jt != (*it).second.end(); ++jt)
    {
      if (jt->first > 100) IsReal.push_back (false);
      else IsReal.push_back (true);
    }
    return IsReal;
  }

  //--------------prendere per ogni cella i punti sopra o sotto la faglia-------------------

  void CProp::getCellPoints (Intersect::GridIntersections_Const_Iterator_Type& it, std::vector <Point3D>& punti, int  ispositive, Point3D const & puntosurf) const
  {
    CPcell cella (M_gridpointer->cell ( (*it).second.i(), (*it).second.j(), (*it).second.k() ) );
    //    int puntis(punti.size());
    for (int ii = 1; ii <= 8; ++ii)
    {
      //devo capire da che parte sta il nodo
      //uso la normale
      Point3D diff = cella.getVertex (ii) - puntosurf;
      if (diff.dot (M_faultpointer->normal (0.5, 0.5) ) *ispositive > 0  ) //sono dalla parte giusta
      {
        punti.push_back (cella.getVertex (ii) );
      }

    }


  }

  Real CProp::setIntArea (Hull const & calimero, gmm::size_type npunti_faglia)const
  {
    Real Area (0);
    for (gmm::size_type i = 0; i < calimero.getNtetra(); ++i)
    {
      std::vector<gmm::size_type> punti_tetra (calimero.getPointsSimplex (i) );
      gmm::size_type cont (0);
      std::vector<gmm::size_type> quali;
      for (gmm::size_type j = 0; j < punti_tetra.size(); ++j)
      {
        if (punti_tetra[j] < npunti_faglia)
        {
          cont += 1;
          quali.push_back (punti_tetra[j]);
        }
      }

      if (cont == 3)
      {
        Triangle t (calimero.getPoint (quali[0]), calimero.getPoint (quali[1]), calimero.getPoint (quali[2]) );
        Area = Area + t.area();
      }
    }

    return Area;
  }

  Point3D CProp::setCG (std::vector<Point3D> const & punti)const
  {
    Point3D CG (0, 0, 0);

    for (gmm::size_type i = 0; i < punti.size(); ++i)
    {
      CG = CG + punti[i];

    }

    CG.x = (1. / Real (punti.size() ) ) * CG.x;
    CG.y = (1. / Real (punti.size() ) ) * CG.y;
    CG.z = (1. / Real (punti.size() ) ) * CG.z;

    return CG;
  }

  Point3D CProp::setCG (Hull const & calimero, gmm::size_type npunti_faglia)const
  {
    Point3D CG (0, 0, 0);
    Real Area (0);
    for (gmm::size_type i = 0; i < calimero.getNtetra(); ++i)
    {
      std::vector<gmm::size_type> punti_tetra (calimero.getPointsSimplex (i) );
      gmm::size_type cont (0);
      std::vector<gmm::size_type> quali;
      for (gmm::size_type j = 0; j < punti_tetra.size(); ++j)
      {
        if (punti_tetra[j] < npunti_faglia)
        {
          cont += 1;
          quali.push_back (punti_tetra[j]);
        }
      }

      if (cont == 3)
      {
        Triangle t (calimero.getPoint (quali[0]), calimero.getPoint (quali[1]), calimero.getPoint (quali[2]) );
        Area = Area + t.area();
        CG = CG + (calimero.getPoint (quali[0]) + calimero.getPoint (quali[1]) + calimero.getPoint (quali[2]) ) * t.area() / 3.;
      }
    }
    return CG / Area;
  }

  Real CProp::setIntd (Hull& calimero, Point3D Aa)
  {
    Real Intd (0);
    Tetra t (calimero.getTetra (0) );
    Point3D N (M_faultpointer->normal (1, 1) );
    std::vector<Real> pesi (t.getGaussWeights() );
    for (gmm::size_type i = 0; i < calimero.getNtetra(); ++i)
    {
      Real Intdprovv (0);
      t = calimero.getTetra (i);
      std::vector<Point3D> nodi (t.getGaussNodes() );
      for (gmm::size_type j = 0; j < nodi.size(); ++j)
      {
        Point3D diff;
        diff = Aa - nodi[j];
        Intdprovv = Intdprovv + pesi[j] * gmm::abs (diff.dot (N) );
      }
      Intd = Intd + Intdprovv * t.volume();
    }
    return Intd;
  }

  std::vector<Point3D> CProp::addPoints4area (std::vector<Point3D> const & puntiarea, std::vector<bool> const & puntiIsReal, Intersect::GridIntersections_Const_Iterator_Type& it)const
  {

    CPcell cella (M_gridpointer->cell ( (*it).second.i(), (*it).second.j(), (*it).second.k() ) );
    gmm::size_type npunti = puntiarea.size();
    std::vector<Point3D> puntiareaNew;
    std::vector<Point3D>::iterator location;
    // Punti che definiscono la superficie bilineare
    // nello spazio fisico
    Point3D A (M_faultpointer->A() );
    Point3D B (M_faultpointer->B() );
    Point3D C (M_faultpointer->C() );
    Point3D D (M_faultpointer->D() );

    // Segmenti che definiscono il bordi del patch nello
    // spazio fisico
    Segment uno (A, B);
    Segment due (B, C);
    Segment tre (C, D);
    Segment quattro (D, A);

    // Loops on all interection points
    for (gmm::size_type i = 0; i < npunti; ++i)
    {
    // Loops on all interection points
      for (gmm::size_type j = i + 1; j < npunti; ++j)
      {
	// If at least one point is vistual/real and we are not 
	// taking the same point
	// @note Point3D::operator == does not comply with the semantic
	// of a equality operator. Is what we want here? 
	// We really want not to consider points that are too nearby? Or we want
	// strict equality?
        if (puntiIsReal[i] != puntiIsReal[j] && ! (puntiarea[i] == puntiarea[j]) )
        {
          Point3D medio;
	  // !Add new points
          Segment linea (puntiarea[i], puntiarea[j]);
          if (linea.intersectTheSegment (uno, medio) )
          {
            puntiareaNew.push_back (medio);
          }
          if (linea.intersectTheSegment (due, medio) )
          {
            puntiareaNew.push_back (medio);
          }
          if (linea.intersectTheSegment (tre, medio) )
          {
            puntiareaNew.push_back (medio);
          }
          if (linea.intersectTheSegment (quattro, medio) )
          {
            puntiareaNew.push_back (medio);
          }
        }
      }
    }
    // Add extrema if in cell
    if ( cella.isIn (A) )
    {
      puntiareaNew.push_back (A);
    }
    if (cella.isIn (B) )
    {
      puntiareaNew.push_back (B);
    }
    if (cella.isIn (C) )
    {
      puntiareaNew.push_back (C);
    }
    if (cella.isIn (D) )
    {
      puntiareaNew.push_back (D);
    }
    
    // a points at the intersection between fracture boundary and 
    // cell faces
    for (gmm::size_type ff = 0; ff < 6; ++ff)
    {
      Point3D medio;

      if (cella.intersectTheFace (uno, ff,   medio) )
      {
        puntiareaNew.push_back (medio);
      }
      if (cella.intersectTheFace (due, ff,   medio) )
      {
        puntiareaNew.push_back (medio);
      }
      if (cella.intersectTheFace (tre, ff,   medio) )
      {
        puntiareaNew.push_back (medio);
      }
      if (cella.intersectTheFace (quattro, ff,   medio) )
      {
        puntiareaNew.push_back (medio);
      }
    }

    // Add all the real points
    for (gmm::size_type i = 0; i < npunti; ++i)
    {
      if (puntiIsReal[i] == true)
      {
        puntiareaNew.push_back (puntiarea[i]);
      }
    }

    return puntiareaNew;
  }

  Real CProp::setIntdist_linea (Hull& calimero,  gmm::size_type npunti_faglia, Fracture::IntFrac intersezione)
  {
    Real Intd (0);
    Real Area (0);
    for (gmm::size_type i = 0; i < calimero.getNtetra(); ++i)
    {
      std::vector<gmm::size_type> punti_tetra (calimero.getPointsSimplex (i) );
      gmm::size_type cont (0);
      std::vector<gmm::size_type> quali;
      for (gmm::size_type j = 0; j < punti_tetra.size(); ++j)
      {
        if (punti_tetra[j] < npunti_faglia)
        {
          cont += 1;
          quali.push_back (punti_tetra[j]);
        }
      }

      if (cont == 3)
      {

        Triangle t (calimero.getPoint (quali[0]), calimero.getPoint (quali[1]), calimero.getPoint (quali[2]) );
        std::vector<Real> pesi (t.getGaussWeights (4) );
        Real Intdprovv (0);
        std::vector<Point3D> nodi (t.getGaussNodes (4) );
        for (gmm::size_type j = 0; j < nodi.size(); ++j)
        {
          Point3D diff;
          diff = intersezione.SMax.A() - nodi[j];

          Intdprovv = Intdprovv + pesi[j] * gmm::abs (diff.dot (intersezione.Normale) ) / intersezione.Normale.norm();

        }
        Area = Area + t.area();
        Intd = Intd + Intdprovv * t.area();
      }
    }

    return Intd / Area;
  }

  Real CProp::setIntdist_linea (std::vector<Point3D> puntiarea, Point3D normale, Fracture::IntFrac intersezione, Intersect::GridIntersections_Const_Iterator_Type& it, bool isPartial)
  {
    Real Intd (0);
    Real Area (0);
    Point3D normale_linea;
    normale_linea = normale.cross (intersezione.SMax.A() - intersezione.SMax.B() );
    CPcell cella (M_gridpointer->cell ( (*it).second.i(), (*it).second.j(), (*it).second.k() ) );
    std::vector<Point3D> puntilinea, puntilineaprovv;

    //seleziono i punti da una parte
    std::vector<Point3D> puntisin;
    for (gmm::size_type i = 0; i < puntiarea.size(); ++i)
    {
      Point3D diff (puntiarea[i] - intersezione.SMax.A() );
      if (diff.dot (normale_linea) > 0)
      {
        puntisin.push_back (puntiarea[i]);
      }

    }
    //seleziono i punti dall'altra parte
    std::vector<Point3D> puntidx;
    for (gmm::size_type i = 0; i < puntiarea.size(); ++i)
    {
      Point3D diff (puntiarea[i] - intersezione.SMax.A() );
      if (diff.dot (normale_linea) <= 0)
      {
        puntidx.push_back (puntiarea[i]);
      }

    }

    Segment SSMax;
    SSMax.setA (20 * intersezione.SMax.A() - 19 * intersezione.SMax.B() );
    SSMax.setB (20 * intersezione.SMax.B() - 19 * intersezione.SMax.A() );

    for (gmm::size_type i = 0; i < puntisin.size(); ++i)
    {
      for (gmm::size_type j = 0; j < puntidx.size(); ++j)
      {

        Point3D A (puntisin[i]), B (puntidx[j]), medio;
        Segment SS (A, B);
        SSMax.intersectTheSegment (SS, medio);
        if (SS.isIn (medio) && cella.isIn (medio) )
        {
          puntilineaprovv.push_back (medio);
        }
      }
    }

    if (puntilineaprovv.size() >= 2)
    {
      Segment SSS (maxSegment (puntilineaprovv) );
      puntilinea.push_back (SSS.A() );
      puntilinea.push_back (SSS.B() );
    }

    if (puntisin.size() > 0)
    {
      for (gmm::size_type ii = 0; ii < puntilinea.size(); ++ii)
      {
        puntisin.push_back (puntilinea[ii]);

      }


      puntisin.push_back (puntisin[0] + normale * intersezione.SMax.length() );

      //chiamo qhull
      Hull guscio1 (puntisin);
      gmm::size_type npunti_faglia = puntisin.size() - 1;

      //calcolo tutto a sin

      for (gmm::size_type i = 0; i < guscio1.getNtetra(); ++i)
      {
        std::vector<gmm::size_type> punti_tetra (guscio1.getPointsSimplex (i) );
        gmm::size_type cont (0);
        std::vector<gmm::size_type> quali;
        for (gmm::size_type j = 0; j < punti_tetra.size(); ++j)
        {
          if (punti_tetra[j] < npunti_faglia)
          {
            cont += 1;
            quali.push_back (punti_tetra[j]);
          }
        }

        if (cont == 3)
        {

          Triangle t (guscio1.getPoint (quali[0]), guscio1.getPoint (quali[1]), guscio1.getPoint (quali[2]) );
          std::vector<Real> pesi (t.getGaussWeights (4) );
          Real Intdprovv (0);
          std::vector<Point3D> nodi (t.getGaussNodes (4) );
          for (gmm::size_type j = 0; j < nodi.size(); ++j)
          {
            Point3D diff;
            diff = intersezione.SMax.A() - nodi[j];

            Intdprovv = Intdprovv + pesi[j] * gmm::abs (diff.dot (intersezione.Normale) ) / intersezione.Normale.norm();

          }
          Area = Area + t.area();
          Intd = Intd + Intdprovv * t.area();
        }
      }
    }


    if (puntidx.size() > 0 )
    {
      for (gmm::size_type ii = 0; ii < puntilinea.size(); ++ii)
      {
        puntidx.push_back (puntilinea[ii]);
      }

      puntidx.push_back (puntidx[0] + normale * intersezione.SMax.length() );

      //chiamo qhull
      Hull guscio2 (puntidx);
      gmm::size_type npunti_faglia = puntidx.size() - 1;
      //calcolo tutto a dx
      for (gmm::size_type i = 0; i < guscio2.getNtetra(); ++i)
      {
        std::vector<gmm::size_type> punti_tetra (guscio2.getPointsSimplex (i) );
        gmm::size_type cont (0);
        std::vector<gmm::size_type> quali;
        for (gmm::size_type j = 0; j < punti_tetra.size(); ++j)
        {
          if (punti_tetra[j] < npunti_faglia)
          {
            cont += 1;
            quali.push_back (punti_tetra[j]);
          }
        }

        if (cont == 3)
        {

          Triangle t (guscio2.getPoint (quali[0]), guscio2.getPoint (quali[1]), guscio2.getPoint (quali[2]) );
          std::vector<Real> pesi (t.getGaussWeights (4) );
          Real Intdprovv (0);
          std::vector<Point3D> nodi (t.getGaussNodes (4) );
          for (gmm::size_type j = 0; j < nodi.size(); ++j)
          {
            Point3D diff;
            diff = intersezione.SMax.A() - nodi[j];

            Intdprovv = Intdprovv + pesi[j] * gmm::abs (diff.dot (intersezione.Normale) ) / intersezione.Normale.norm();

          }
          Area = Area + t.area();
          Intd = Intd + Intdprovv * t.area();
        }
      }
    }


    return Intd / Area;

  }

} // namespace Geometry
