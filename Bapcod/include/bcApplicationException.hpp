/**
 *
 * This file bcApplicationException.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCAPPLICATIONEXCEPTION_HPP_
#define BCAPPLICATIONEXCEPTION_HPP_

#include <string>
#include <iostream>

class GlobalException;

/**
 * This exception class must be used by the applications.
 * It provides also all the function call strace.
 */
class ApplicationException
{
  GlobalException* _exception;
public:

  /**
   * Constructor.
   * @param message Exception message.
   * @param printingForced Is the message is displayed or not?
   * @param os OutputStream
   */
  ApplicationException(std::string message="", bool printingForced = true, std::ostream & os = std::cerr);
  virtual ~ApplicationException() throw ();

  const std::string& message() const;

  std::ostream& printTrace(std::ostream& os = std::cout) const;
};

std::ostream& operator<<(std::ostream& os, const ApplicationException& exception);


#endif /* BCAPPLICATIONEXCEPTION_HPP_ */
