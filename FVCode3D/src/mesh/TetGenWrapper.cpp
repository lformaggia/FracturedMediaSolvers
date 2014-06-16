/*!
 *	@file TetGenWrapper.cpp
 *	@brief Class for tetrahedralizing a polyhedron (definitions).
 */

#include "tetgen.h"
#include "geometry/Operations.hpp"
#include "mesh/TetGenWrapper.hpp"

namespace FVCode3D
{

TetGenWrapper::TetGenWrapper(std::vector<Point3D> nodes, std::vector< std::vector<UInt> > faces):
	M_inNodes(nodes), M_faces(faces), M_volume(0.) {}

void TetGenWrapper::generateMesh()
{
	tetgenio in, out;
	tetgenio::facet *f;
	tetgenio::polygon *p;
	Point3D tmp;

	// All indices start from 1.
	in.firstnumber = 0;

	// Set nodes
	in.numberofpoints = M_inNodes.size();
	in.pointlist = new REAL[in.numberofpoints * 3];

	for(Int i=0; i<in.numberofpoints; ++i)
	{
		in.pointlist[3*i]   = M_inNodes[i].x();
		in.pointlist[3*i+1] = M_inNodes[i].y();
		in.pointlist[3*i+2] = M_inNodes[i].z();
	}

	// Set facets
	in.numberoffacets = M_faces.size();
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];

	for(Int i=0; i<in.numberoffacets; ++i)
	{
		f = &in.facetlist[i];
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
		f->holelist = NULL;
		p = &f->polygonlist[0];
		p->numberofvertices = M_faces[i].size();
		p->vertexlist = new int[p->numberofvertices];
		for(Int j=0; j<p->numberofvertices; ++j)
		{
			p->vertexlist[j] = M_faces[i][j];
		}
	}

	// Set 'in.facetmarkerlist'
	for(Int i=0; i<in.numberoffacets; ++i)
		in.facetmarkerlist[i] = 0;

	char option[] = "pYQ";
	tetrahedralize(option, &in, &out);

	// Get the output nodes
	M_outNodes.resize(out.numberofpoints);
	for(Int i=0; i<out.numberofpoints; ++i)
	{
		tmp.setValues(out.pointlist[3*i], out.pointlist[3*i+1], out.pointlist[3*i+2]);
		M_outNodes[i] = tmp;
	}

	// Get the output elements
	M_elements.resize(out.numberoftetrahedra);
	for(Int i=0; i<out.numberoftetrahedra; ++i)
	{
		M_elements[i].resize(4);
		M_elements[i][0] = out.tetrahedronlist[4*i];
		M_elements[i][1] = out.tetrahedronlist[4*i+1];
		M_elements[i][2] = out.tetrahedronlist[4*i+2];
		M_elements[i][3] = out.tetrahedronlist[4*i+3];
	}
}

Real TetGenWrapper::computeVolume()
{
	M_volume = 0.;
	const UInt nEle = M_elements.size();
	std::vector<Point3D> tmp(4);

	for(UInt i=0; i<nEle; ++i)
	{
		tmp[0] = M_outNodes[M_elements[i][0]];
		tmp[1] = M_outNodes[M_elements[i][1]];
		tmp[2] = M_outNodes[M_elements[i][2]];
		tmp[3] = M_outNodes[M_elements[i][3]];
		M_volume += tetrahedronVolume(tmp);
	}

	return M_volume;
}

Point3D TetGenWrapper::computeCenterOfMass()
{
	M_centerOfMass.setValues(0.,0.,0.);
	const UInt nEle = M_elements.size();
	std::vector<Point3D> tmp(4);
	Real tmpVol, totVol = 0.;

	for(UInt i=0; i<nEle; ++i)
	{
		tmp[0] = M_outNodes[M_elements[i][0]];
		tmp[1] = M_outNodes[M_elements[i][1]];
		tmp[2] = M_outNodes[M_elements[i][2]];
		tmp[3] = M_outNodes[M_elements[i][3]];
		tmpVol = tetrahedronVolume(tmp);
		totVol += tmpVol;
		M_centerOfMass += tmpVol * (tmp[0] + tmp[1] + tmp[2] + tmp[3]) / 4.;
	}

	M_centerOfMass /= totVol;

	return M_centerOfMass;
}

const std::vector<Point3D> TetGenWrapper::getElement(const UInt i) const
{
	std::vector<Point3D> points(4);

	points[0] = M_outNodes[M_elements[i][0]];
	points[1] = M_outNodes[M_elements[i][1]];
	points[2] = M_outNodes[M_elements[i][2]];
	points[3] = M_outNodes[M_elements[i][3]];

	return points;
}

}// namespace FVCode3D
