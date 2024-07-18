/**
 *
 * This file bcGlobalException.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcGlobalException.hpp"
//#include <boost/stacktrace.hpp>

GlobalException::GlobalException(std::string message, bool printingForced, std::ostream & os) :
    _message(message)
{
  if (printingForced)
      os << *this;
}

GlobalException::~GlobalException() throw ()
{
}

const std::string& GlobalException::message() const
{
  return _message;
}

std::ostream & operator <<(std::ostream & os, const GlobalException & exception)
{
  os << exception.message() << std::endl;
  //os << bcoost::stacktrace::stacktrace();
  return os;
}

