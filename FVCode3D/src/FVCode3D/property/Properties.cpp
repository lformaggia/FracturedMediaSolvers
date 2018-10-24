/*!
 *  @file Properties.cpp
 *  @brief Classes that handle the properties (definitions).
 */

#include <FVCode3D/property/Properties.hpp>
#include <FVCode3D/property/Permeability.hpp>
#include <FVCode3D/property/PermeabilityFactory.hpp>
#include <FVCode3D/mesh/Mesh3D.hpp>
#include <FVCode3D/utility/readPermeabilities.hpp>
#include <FVCode3D/utility/readOtherData.hpp>
#include <iostream>

namespace FVCode3D
{

void PropertiesMap::setPropertiesOnMatrix(const Mesh3D & mesh, const Real porosity, const PermPtr_Type & permeability) // @suppress("Member declaration not found")
{
    std::set<UInt> zones;
    for(std::map<UInt,Mesh3D::Cell3D>::const_iterator it = mesh.getCellsMap().begin(); it != mesh.getCellsMap().end(); ++it)
        if(it->second.getZoneCode() > 0)
            zones.insert(it->second.getZoneCode());

    for(std::set<UInt>::const_iterator it = zones.begin(); it != zones.end(); ++it)
    {
        getProperties(*it).M_aperture = 1.;
        getProperties(*it).M_porosity = porosity;
        getProperties(*it).M_permeability = permeability;
    }
}

void PropertiesMap::setPropertiesOnFractures(const Mesh3D & mesh, const Real aperture, const Real porosity, const PermPtr_Type & permeability) // @suppress("Member declaration not found")
{
    std::set<UInt> zones;
    for(std::map<UInt,Mesh3D::Facet3D>::const_iterator it = mesh.getFacetsMap().begin(); it != mesh.getFacetsMap().end(); ++it)
        if(it->second.getZoneCode() > 0)
            zones.insert(it->second.getZoneCode());

    for(std::set<UInt>::const_iterator it = zones.begin(); it != zones.end(); ++it)
    {
        getProperties(*it).M_aperture = aperture;
        getProperties(*it).M_porosity = porosity;
        getProperties(*it).M_permeability = permeability;
    }
}

void PropertiesMap::setPropertiesOnFractures(const Mesh3D & mesh, // @suppress("Member declaration not found")
        const Utility::fractureDataList & porosityList,
        Utility::FracturePermeabilityData const & permeabilities)
{
    static bool onlyone(false);
    std::set<UInt> zones;
    for(std::map<UInt,Mesh3D::Facet3D>::const_iterator it = mesh.getFacetsMap().begin();
            it != mesh.getFacetsMap().end(); ++it)
        if(it->second.getZoneCode() > 0)
           zones.insert(it->second.getZoneCode());
    // I need default permeability!
    auto defaultPermeability=Utility::getPermeability(permeabilities,0);
    // I need default porosity.
    auto defaultFractureData = Utility::getFractureData(porosityList,0);
    PermeabilityFactory & permeabilityFactory= PermeabilityFactory::Instance();

    for(std::set<UInt>::const_iterator it = zones.begin(); it != zones.end(); ++it)
    {
        auto perm=defaultPermeability;
        auto poroAndAperture=defaultFractureData;

        try
        {
            perm=Utility::getPermeability(permeabilities,*it);
        }
        catch (...)
        {
            if(! onlyone){
            onlyone=true;
            std::cerr<<"Cannot find permeability associated to fracture zone "<<*it<<std::endl;
            std::cerr<<"Using default"<<std::endl;
            }
            //perm=defaultPermeability;
        }
        try
        {
            poroAndAperture=Utility::getFractureData(porosityList,*it);
        }
        catch(...)
        {
            if(!onlyone)
              {
            std::cerr<<"Cannot find porosity associated to fracture zone "<<*it<<std::endl;
            std::cerr<<"Using default"<<std::endl;
              }
            //poroAndAperture=defaultFractureData;
        }
        // Now build the real stuff

        PermPtr_Type p_ptr;
        auto ptype=std::get<0>(perm);
        auto pvalues=std::get<1>(perm);
        switch (ptype)
        {
            case 0:
                p_ptr = permeabilityFactory.getProduct("ScalarPermeability");
                p_ptr->setPermeability (pvalues[0],0);
                break;
            case 1:
                p_ptr = permeabilityFactory.getProduct("DiagonalPermeability");
                p_ptr->setPermeability (pvalues[0],0,0);
                p_ptr->setPermeability (pvalues[1],1,1);
                p_ptr->setPermeability (pvalues[2],2,2);
                break;
            case 2:
                p_ptr = permeabilityFactory.getProduct("SymTensorPermeability");
                p_ptr->setPermeability (pvalues[0],0,0);
                p_ptr->setPermeability (pvalues[1],0,1);
                p_ptr->setPermeability (pvalues[2],0,2);
                p_ptr->setPermeability (pvalues[3],1,1);
                p_ptr->setPermeability (pvalues[4],1,2);
                p_ptr->setPermeability (pvalues[5],2,2);
                break;
        }

        getProperties(*it).M_aperture = std::get<1>(poroAndAperture);
        getProperties(*it).M_porosity = std::get<0>(poroAndAperture);
        getProperties(*it).M_permeability = p_ptr;
    }
}

void PropertiesMap::setPropertiesOnMatrix(const Mesh3D & mesh, // @suppress("Member declaration not found")
        const Utility::bulkDataList & porosityList,
        Utility::BulkPermeabilityData const & permeabilities)
{
    static bool onlyone(false);
    std::set<UInt> zones;
    for(std::map<UInt,Mesh3D::Cell3D>::const_iterator it = mesh.getCellsMap().begin();
            it != mesh.getCellsMap().end(); ++it)
        if(it->second.getZoneCode() > 0)
           zones.insert(it->second.getZoneCode());
    constexpr Real defaultAperture=1.0;
    // I need default permeability!
    auto defaultPermeability=Utility::getPermeability(permeabilities,0);
    // I need default porosity.
    auto defaultPorosity = Utility::getBulkData(porosityList,0);
    PermeabilityFactory & permeabilityFactory= PermeabilityFactory::Instance();

    for(std::set<UInt>::const_iterator it = zones.begin(); it != zones.end(); ++it)
    {
        Utility::Permeability perm;
        Real poro=defaultPorosity;

        try
        {
            perm=Utility::getPermeability(permeabilities,*it);
        }
        catch (...)
        {
            if (!onlyone)
              {
            std::cerr<<"Cannot find permeability associated to matrix zone "<<*it<<std::endl;
            std::cerr<<"Using default (will be printed only once)"<<std::endl;
            onlyone=true;
              }
            perm=defaultPermeability;
        }
        try
        {
            poro=Utility::getBulkData(porosityList,*it);
        }
        catch(...)
        {
            if (! onlyone){
            std::cerr<<"Cannot find porosity associated to matrix zone "<<*it<<std::endl;
            std::cerr<<"Using default"<<std::endl;
            }
            poro=defaultPorosity;
        }
        // Now build the real stuff

        PermPtr_Type p_ptr;
        auto ptype=std::get<0>(perm);
        auto pvalues=std::get<1>(perm);
        switch (ptype)
        {
            case 0:
                p_ptr = permeabilityFactory.getProduct("ScalarPermeability");
                p_ptr->setPermeability (pvalues[0],0);
                break;
            case 1:
                p_ptr = permeabilityFactory.getProduct("DiagonalPermeability");
                p_ptr->setPermeability (pvalues[0],0,0);
                p_ptr->setPermeability (pvalues[1],1,1);
                p_ptr->setPermeability (pvalues[2],2,2);
                break;
            case 2:
                p_ptr = permeabilityFactory.getProduct("SymTensorPermeability");
                p_ptr->setPermeability (pvalues[0],0,0);
                p_ptr->setPermeability (pvalues[1],0,1);
                p_ptr->setPermeability (pvalues[2],0,2);
                p_ptr->setPermeability (pvalues[3],1,1);
                p_ptr->setPermeability (pvalues[4],1,2);
                p_ptr->setPermeability (pvalues[5],2,2);
                break;
        }

        getProperties(*it).M_aperture = defaultAperture;
        getProperties(*it).M_porosity = poro;
        getProperties(*it).M_permeability = p_ptr;
    }
}

} // namespace FVCode3D
