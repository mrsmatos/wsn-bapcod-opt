/**
 *
 * This file bcParameterParserC.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCPARAMETERPARSERC_H_
#define BCPARAMETERPARSERC_H_

#include <boost/program_options.hpp>
#include <boost/any.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>

#include "bcModelParameterC.hpp"

/**
 * This class permits to parse a file containing some parameters.
 */
class ParameterParser
{
protected:
  boost::program_options::options_description _configFileOptions;
  boost::program_options::options_description _cmdLineOptions;
  boost::program_options::options_description _visibleOptions;

  std::string _parameterFileName;
  bool _isAlreadyParsed;

public:

  /**
   * Constructor
   * @param parameterFileName filename of your parameter.
   * If parameterFileName is not passed, the parameter searching will in the same file as BaPCod parameters.
   */
  ParameterParser(std::string parameterFileName = std::string(""));

  virtual ~ParameterParser();

  /**
   * Parse the parameter file
   * @param argc nb arguments
   * @param argv arguments
   * @throw GlobalException if the file is already parsed.
   */
  virtual void parse(int argc, char** argv);

  virtual void parseFile(std::string parameterFileName);

  void setParameterFileName(const std::string & fileName);

  std::string getParameterFileName() const;

  /**
   * Add an application parameter (an application parameter for example) to the
   * program parameter.
   *
   * @param name parameter name (inside the configuration file if
   *        configOption)
   * @param value pointer of the value to the parameter
   * @param defaultValue default value of the parameter
   * @param description description of the parameter
   * @param configOption is it a configuration parameter (inside the file)?
   * @param hideOption is option displaying (when we call the application with
   *        "help" parameter)?
   */
  template<typename T>
  void addApplicationParameter(std::string name, T* value,
                               const T& defaultValue,
                               std::string description,
                               bool configOption = true,
                               bool hideOption = false)
  {
    _cmdLineOptions.add_options()(name.c_str(),boost::program_options::value<T>(value)->default_value(defaultValue),
                                  description.c_str());

    if (configOption || hideOption)
      {
        _configFileOptions.add_options()(name.c_str(),
                                         boost::program_options::value<T>(value)->default_value(defaultValue),
                                         description.c_str());
      }

    if (!hideOption)
      {
        _visibleOptions.add_options()(name.c_str(),
                                      boost::program_options::value<T>(value)->default_value(defaultValue),
                                      description.c_str());
      }
  }

  /**
   * Add an application parameter (an application parameter for example) to the
   * program parameter.
   *
   * @param name parameter name (inside the configuration file if
   *        configOption)
   * @param applicationParameter: parameter application
   * @param configOption is it a configuration parameter (inside the file)?
   * @param hideOption is option displaying (when we call the application with
   *        "help" parameter)?
   */
  template<typename T, typename U>
    void addApplicationParameter(ApplicationParameter<T, U>& applicationParameter,
                                 bool configOption = true,
                                 bool hideOption = false)
    {
      addApplicationParameter(applicationParameter.name(),
                              applicationParameter.valuePtr(), applicationParameter.defaultValue(),
                              applicationParameter.description(), configOption, hideOption);
    }

private:

  void getBapcodParameterFilename(int argc, char** argv);
};


#endif /* BCPARAMETERPARSERC_H_ */
