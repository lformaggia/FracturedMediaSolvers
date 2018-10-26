#include <FVCode3D/utility/functionLoader.hpp>
#include <stdexcept>
#include <unordered_map>
#include <FVCode3D/core/functions.hpp>
FVCode3D::FunctionLoader::FunctionLoader ()
{
  // Insert predefined functions
  this->M_cached.insert({"fZero",FVCode3D::fZero});
  this->M_cached.insert({"fOne", FVCode3D::fOne} );
  this->M_cached.insert({"fTwo", FVCode3D::fTwo} );
}

FVCode3D::FunctionLoader::FunctionLoader (std::string fileName): FunctionLoader()
{
  this->load(fileName);
}


void FVCode3D::FunctionLoader::load(std::string fileName)
{
  this->release();
  auto lib = dlopen(fileName.c_str(),RTLD_NOW);
  if (!lib)
    {
    lib=nullptr;
    throw std::runtime_error(std::string("Dynamic library named ")+fileName+std::string(" not found!"));
    }
  this->M_lib_handle=lib;
}


FVCode3D::Func
FVCode3D::FunctionLoader::getFunction (std::string functionName) const
{
  // see if present in cache
  auto search = this->M_cached.find(functionName);
  if (search != this->M_cached.end())
    {
      return search->second;
    }
  else
    {
      if (this->M_lib_handle == nullptr)
        throw std::runtime_error(std::string("Dynamic library for function ")+functionName+std::string(" not loaded!"));
      auto funp = reinterpret_cast<FVCode3D::Func *>(dlsym(this->M_lib_handle,functionName.c_str()));
      char * error = dlerror();
      if (error != nullptr || funp == nullptr)
        throw std::runtime_error(std::string("Function ")+functionName+std::string(" not found.\n Error=")+std::string(error));
      // insert to simplify search if needed next time
      this->M_cached.insert(std::make_pair(functionName,*funp));
      return FVCode3D::Func(*funp);
    }
}
