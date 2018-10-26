/*
 * functionLoader.hpp
 *
 *  Created on: Oct 25, 2018
 *      Author: forma
 */

#ifndef SUBPROJECTS__FVCODE3D_SRC_FVCODE3D_UTILITY_FUNCTIONLOADER_HPP_
#define SUBPROJECTS__FVCODE3D_SRC_FVCODE3D_UTILITY_FUNCTIONLOADER_HPP_
#include <FVCode3D/core/TypeDefinition.hpp>
#include <dlfcn.h>
#include <memory>
#include <string>
#include <unordered_map>
namespace FVCode3D
{
  class FunctionLoader
  {
  public:
    FunctionLoader();
    //! Creates the loader and opens the library
    /*!
     * @par fileName a valid shared library containing Func objects
     * @throw runtime_error if library not present
     */
    FunctionLoader(std::string fileName);
    //! Loads the library
    /*!
     * @throw runtime_error if library not present
     */
    void load(std::string fileName);
    // releases the library
    void release()
    {
      if (M_lib_handle != nullptr)
        {
          dlclose(M_lib_handle);
          M_lib_handle=nullptr;
        }
    }
    //! extract function
    /*!
     * @par functionName the name of the function
     * @return a Func object with the function
     * @throw runtime__error runtime exception if name not found
     */
    Func getFunction(std::string functionName) const;
  private:
    void * M_lib_handle{nullptr};
    mutable std::unordered_map<std::string,FVCode3D::Func> M_cached;
  };
}


#endif /* SUBPROJECTS__FVCODE3D_SRC_FVCODE3D_UTILITY_FUNCTIONLOADER_HPP_ */
