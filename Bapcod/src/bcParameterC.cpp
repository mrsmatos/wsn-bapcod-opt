/**
 *
 * This file bcParameterC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

// /**
//  *  BaPCod  -  a generic Branch-And-Price Code
//  *
//  *  @file   bcParameterC.cpp
//  *
//  *  @author Francois VANDERBECK <Francois.Vanderbeck@math.u-bordeaux1.fr>,
//  *          Romain LEGUAY <romain.leguay@math.u-bordeaux1.fr>
//  *  All Rights Reserved. [LICENCE]
//  *
//  *  @brief Parameters read by boost program_option library.
//  *
//  *   Created on: 23 avr. 2012
//  */

#include "bcParameterC.hpp"
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>


StrongBranchingPhaseParameter::StrongBranchingPhaseParameter(const std::string & str):
  _active(false), _exact(false), _maxNumOfCandidates(0), _maxNumOfColGenIterations(0),
  _minLevelOfSpRestriction(0), _minNumCutRounds(0), _maxNumCutRounds(10000),
  _doRedCostFixingAndEnumeration(true), _logPrintFrequency(0),  _treeSizeRatioToStop(1.0)
{
}

StrongBranchingPhaseParameter::~StrongBranchingPhaseParameter()
{
}

std::ostream & StrongBranchingPhaseParameter::print(std::ostream & os) const
{
  if (!_active)
    os << "not active";
  else
    {
      os << " max#cand. = " << _maxNumOfCandidates;
      os << " max#cg.iters = " << _maxNumOfColGenIterations;
      os << " min.lvl.sp.restr. = " << _minLevelOfSpRestriction;
      os << " min#cut.rounds = " << _minNumCutRounds;
      os << " max#cut.rounds = " << _maxNumCutRounds;
      os << " red.cost.fix&enum. = " << _doRedCostFixingAndEnumeration;
      os << " tree.size.ratio = " << _treeSizeRatioToStop;
    }
  return os;
}

void StrongBranchingPhaseParameter::validate(boost::any& v, const std::vector<std::string>& values,
                                             std::string s)
{
  boost::char_separator<char> sep(s.c_str());
  boost::tokenizer<boost::char_separator<char> > tok(values[0], sep);

  /// default values
  _active = false;
  _maxNumOfCandidates = 1;
  _treeSizeRatioToStop = 1.0;
  _maxNumOfColGenIterations = 10000000;
  _minLevelOfSpRestriction = 0;
  _minNumCutRounds = 0;
  _maxNumCutRounds = 1000000;
  _doRedCostFixingAndEnumeration = true;

  /// Validate and insert each string from the command line
  int counter = 0;
#if (BOOST_VERSION > 104000)
  BOOST_FOREACH(std::string element, tok)
  {
#else
  for(boost::tokenizer<boost::char_separator<char> > ::iterator i = tok.begin(); i != tok.end(); i++)
    {
      std::string element = *i;
#endif //(BOOST_VERSION > 104000)
    try
      {
        switch (counter++) {
          case 0:
            std::istringstream(element) >> std::boolalpha >> _exact;
            break;

          case 1:
            _maxNumOfCandidates = std::atoi(element.c_str());
            break;

          case 2:
            _treeSizeRatioToStop = std::atof(element.c_str());
            break;

          case 3:
            _maxNumOfColGenIterations = std::atoi(element.c_str());
            break;

          case 4:
            _minLevelOfSpRestriction = std::atoi(element.c_str());
            break;

          case 5:
            _minNumCutRounds = std::atoi(element.c_str());
            break;

          case 6:
            _maxNumCutRounds = std::atoi(element.c_str());
            break;

          case 7:
            std::istringstream(element) >> std::boolalpha >> _doRedCostFixingAndEnumeration;
            break;

          case 8:
            _logPrintFrequency = std::atoi(element.c_str());
            break;

          default:
            break;
          }
      }
    catch (const std::exception& e)
      {
        throw boost::program_options::invalid_option_value(
                e.what()
                + std::string(
                  ": Invalid value." + element + std::string("\n")));
      }
    _active = true;
  }

  ///   Use the CRTP to cast, and then clone the container in a boost::any
  v = boost::any(*this);
}

void StrongBranchingPhaseParameter::setExact()
{
    _active = true;
    _exact = true;
    _maxNumOfCandidates = 1;
    _maxNumOfColGenIterations = 10000;
    _minLevelOfSpRestriction = 0;
    _minNumCutRounds = 0;
    _maxNumCutRounds = 10000;
    _doRedCostFixingAndEnumeration = true;
    _logPrintFrequency = 1;
    _treeSizeRatioToStop = 1.0;
}

void StrongBranchingPhaseParameter::setNonActive()
{
    _active = false;
}

void StrongBranchingPhaseParameter::setNonExact(const int maxNumOfCandidates, const int maxNumOfColGenIterations,
                                                const int minLevelOfSpRestriction, const int minNumCutRounds,
                                                const int maxNumCutRounds, const bool doRedCostFixingAndEnumeration,
                                                const int logPrintFrequency, const double treeSizeRatioToStop)
{
    _active = true;
    _exact = false;
    _maxNumOfCandidates = maxNumOfCandidates;
    _maxNumOfColGenIterations = maxNumOfColGenIterations;
    _minLevelOfSpRestriction = minLevelOfSpRestriction;
    _minNumCutRounds = minNumCutRounds;
    _maxNumCutRounds = maxNumCutRounds;
    _doRedCostFixingAndEnumeration = doRedCostFixingAndEnumeration;
    _logPrintFrequency = logPrintFrequency;
    _treeSizeRatioToStop = treeSizeRatioToStop;
}

