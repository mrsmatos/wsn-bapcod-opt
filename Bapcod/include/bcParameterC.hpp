/**
 *
 * This file bcParameterC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCPARAMETERC_H_
#define BCPARAMETERC_H_

#include <boost/program_options.hpp>
#include <boost/any.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/static_assert.hpp>

#if(BOOST_VERSION > 104000)
#include <boost/foreach.hpp>
#endif //(BOOST_VERSION > 104000)

#include <iostream>
#include <fstream>
#include <iterator>

/// Special CRTP that control parameter of class T
template<class T, typename U = int>
class Parameter
{
public:
  virtual void validate(boost::any & v,
      const std::vector<std::string> & values, std::string r = "\\d+")
  {
    boost::regex regexp(r);

    /// Make sure no previous assignment to 'values' was made.
    boost::program_options::validators::check_first_occurrence(v);

    /**
     * Extract the first string from 'values'. If there is more than
     * one string, it's an error, and exception will be thrown.
     */
    const std::string& s =
        boost::program_options::validators::get_single_string(values);

    /**
     * Do regex match and convert the interesting part to
     * int.
     */
    boost::smatch match;
    if (boost::regex_match(s, match, regexp))
      {
        ///     Use the CRTP to cast, and then clone the derived object in a boost::any
        v = boost::any(T(boost::lexical_cast<U>(match.str())));
        return;
      }
    else
      {
        throw boost::program_options::invalid_option_value(
            "Invalid value: " + s);
      }
  }

  virtual ~Parameter() {}
};



/**
 * This class allows to overload custom validator
 * for std::vector.
 * Actualy the multitoken mechanism is empty with config
 * file for std::vector, so that's the only way to overload
 * the validate function
 * Each token in the config file must be separated by a space.
 * Be carefull with the template argument
 *
 * template T is the type of the elements in the container, int, char ...
 * template C the type of the container, set, vector ...
 * template U must be the derived class itself for the CRTP.
 */
template<typename T, typename C, typename U>
class MultitokenParameter : public C
{
public:
  virtual ~MultitokenParameter() {}
  // Function which validates additional tokens from command line.
  virtual void validate(boost::any& v, const std::vector<std::string>& values,
      std::string s = " ")
  {
    boost::char_separator<char> sep(s.c_str());
    boost::tokenizer<boost::char_separator<char> > tok(values[0], sep);

    /// Validate and insert each string from the command line
#if(BOOST_VERSION > 104000)
    BOOST_FOREACH(std::string element, tok)
    {
#else
      for(boost::tokenizer<boost::char_separator<char> > ::iterator i = tok.begin(); i != tok.end(); i++)
        {
          std::string element = *i;
#endif //(BOOST_VERSION > 104000)
          try
          {
              validate_one(element);
          }
          catch (const std::exception& e)
          {
              throw boost::program_options::invalid_option_value(
                  e.what()
                  + std::string(
                      ": Invalid value." + element + std::string("\n")));
          }
        }

      ///   Use the CRTP to cast, and then clone the container in a boost::any
      v = *(dynamic_cast<U*>(this));
    }

private:
    /// This method add an element if we can
    virtual void
    validate_one(std::string element) = 0;
};

  /**
   * This class allow the std::vector to work with multitoken arguments
   * with boost::program_option both at command line and in config file.
   * So we define the validate_one method of MultitokenParameter.
   */
  template<typename T>
  class VectorParameter : public MultitokenParameter<T, std::vector<T>,
  VectorParameter<T> >
  {
  public:
    virtual ~VectorParameter() {}
    virtual void validate_one(std::string token)
    {
      T element = static_cast<T>(token);
      std::vector<T>::push_back(element);
    }
  };

  /**
   * This class allow the std::vector to work with multitoken arguments
   * with boost::program_option both at command line and in config file.
   * So we define the validate_one method of MultitokenParameter.
   */
  template<typename T>
  class SetParameter : public MultitokenParameter<T, std::set<T>,
  SetParameter<T> >
  {
  public:
    virtual ~SetParameter() {}
    virtual void validate_one(std::string element)
    {
      std::set<T>::insert(boost::lexical_cast<T>(element));
    }
  };

class StrongBranchingPhaseParameter
{
  bool _active;
  bool _exact;
  int _maxNumOfCandidates;
  int _maxNumOfColGenIterations;
  int _minLevelOfSpRestriction;
  int _minNumCutRounds;
  int _maxNumCutRounds;
  bool _doRedCostFixingAndEnumeration;
  int _logPrintFrequency;
  double _treeSizeRatioToStop;
  
public:

  StrongBranchingPhaseParameter(const std::string & str = "");

  virtual ~StrongBranchingPhaseParameter();

  std::ostream & print(std::ostream & os) const;

  virtual void validate(boost::any& v, const std::vector<std::string>& values,
                        std::string s = " ");

  inline bool active() const { return _active; }
  inline bool exact() const { return _exact; }
  inline int maxNumOfCandidates() const { return _maxNumOfCandidates; }
  inline int maxNumOfColGenIterations() const { return _maxNumOfColGenIterations; }
  inline int minLevelOfSpRestriction() const { return _minLevelOfSpRestriction; }
  inline int minNumCutRounds() const { return _minNumCutRounds; }
  inline int maxNumCutRounds() const { return _maxNumCutRounds; }
  inline bool doRedCostFixingAndEnumeration() const { return _doRedCostFixingAndEnumeration; }
  inline int logPrintFrequency() const { return _logPrintFrequency; }
  inline double treeSizeRatioToStop() const { return _treeSizeRatioToStop; }

  void setExact();
  void setNonActive();
  void setNonExact(const int maxNumOfCandidates, const int maxNumOfColGenIterations, const int minLevelOfSpRestriction,
                   const int minNumCutRounds, const int maxNumCutRounds, const bool doRedCostFixingAndEnumeration,
                   const int logPrintFrequency, const double treeSizeRatioToStop);
};

inline std::ostream& operator<<(std::ostream& os, const StrongBranchingPhaseParameter & that)
{
  return that.print(os);
}

#endif /* BCPARAMETERC_H_ */
