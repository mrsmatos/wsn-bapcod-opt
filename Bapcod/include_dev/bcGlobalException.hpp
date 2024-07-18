/**
 *
 * This file bcGlobalException.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCGLOBALEXCEPTION_H_
#define BCGLOBALEXCEPTION_H_

#include <string>
#include <iostream>

/**
 * It's a global exception class.
 * The main purpose of this class is to be inherit by other specific exception
 * classes.
 */
class GlobalException
{
protected:
  std::string _message;

public:
  GlobalException(std::string message="", bool printingForced = true, std::ostream & os = std::cerr);
  virtual ~GlobalException() throw ();

  const std::string& message() const;
};

std::ostream& operator<<(std::ostream& os, const GlobalException& exception);

#endif /* BCGLOBALEXCEPTION_H_ */
