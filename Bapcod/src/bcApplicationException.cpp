/**
 *
 * This file bcApplicationException.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcApplicationException.hpp"
#include "bcGlobalException.hpp"

ApplicationException::ApplicationException(std::string message, bool printingForced, std::ostream & os)
{
  _exception = new GlobalException(message, printingForced, os);
  //if the printing force is false, we still print the message.
  if(!printingForced)
    {
      os << message << std::endl;
    }
}

ApplicationException::~ApplicationException() throw ()
{
  delete _exception;
}

const std::string& ApplicationException::message() const
{
  return _exception->message();
}

std::ostream& ApplicationException::printTrace(std::ostream& os) const
{
  //os << boost::stacktrace::stacktrace();
  return os;
}


std::ostream& operator<<(std::ostream& os, const ApplicationException& exception)
{
  os << exception.message() << std::endl;
  exception.printTrace(os);
  return os;
}
