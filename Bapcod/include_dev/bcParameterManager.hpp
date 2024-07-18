/**
 *
 * This file bcParameterManager.hpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#ifndef BCPARAMETERMANAGER_HPP_
#define BCPARAMETERMANAGER_HPP_
#include <boost/program_options.hpp>

#include "bcModelParameterC.hpp"

class ParameterManager
{
private:
  boost::program_options::options_description _configFileOptions;
  boost::program_options::options_description _cmdLineOptions;
  boost::program_options::options_description _visibleOptions;
  boost::program_options::options_description _genericOptions;

public:
  ParameterManager();
  virtual ~ParameterManager();

  /**
   * Parse the parameter file
   * @param argc nb arguments
   * @param argv arguments
   */
  virtual void parse(int argc, char** argv, std::string& filename);

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
  void addParameter(const std::string& name, T* value,
      const T& defaultValue,
      const std::string& description,
      bool configOption = true,
      bool hideOption = true,
      bool genericOption = false)
  {
    if (genericOption)
      {
        addParameterTo(_genericOptions, name, value, defaultValue, description);
      }
    else
      {
        addParameterTo(_cmdLineOptions, name, value, defaultValue, description);
        if (configOption || hideOption)
          {
            addParameterTo(_configFileOptions, name, value, defaultValue,
                description);
          }

        if (!hideOption)
          {
            addParameterTo(_visibleOptions, name, value, defaultValue,
                description);
          }
      }
  }

  /**
   * Add an application parameter (an application parameter for example) to the
   * program parameter.
   *
   * @param applicationParameter: parameter application
   * @param configOption: is it a configuration parameter (inside the file)?
   * @param hideOption: is option displaying (when we call the application with
   *        "help" parameter)?
   * @param genericOption: is it a generic Option (like --help, --version...)?
   */
  template<typename T, typename U>
  void addParameter(ApplicationParameter<T, U>& applicationParameter,
                    bool configOption = true,
                    bool hideOption = true,
                    bool genericOption = false)
  {
    addParameter(applicationParameter.name(),
        applicationParameter.valuePtr(), applicationParameter.defaultValue(),
        applicationParameter.description(), configOption, hideOption,
        genericOption);
  }

  /**
   * Add an advanced parameter (can be a multitoken or not).
   *
   * @param name: parameter's name (inside the configuration file if
   *        configOption)
   * @param value: parameter's value pointer
   * @param defaultValue: parameter's default value
   * @param description: parameter's description.
   * @param isMultitoken: does the parameter contain multiple values?
   * @param configOption: is it a configuration parameter (inside the file)?
   * @param hideOption: is option displaying (when we call the application with
   *        "help" parameter)?
   * @param genericOption: is it a generic Option (like --help, --version...)?
   */
  template<typename T>
  void addAdvancedParameter(const std::string& name,
                            T* value,
                            const std::string& description,
                            bool isMultitoken,
                            bool configOption = true,
                            bool hideOption = true,
                            bool genericOption = false)
  {
    if (genericOption)
      {
        addParameterTo(_genericOptions, name, value,
            description, isMultitoken);
      }
    else
      {
        addParameterTo(_cmdLineOptions, name, value,
            description, isMultitoken);
        if (configOption || hideOption)
          {
            addParameterTo(_configFileOptions, name, value,
                description, isMultitoken);
          }

        if (!hideOption)
          {
            addParameterTo(_visibleOptions, name, value,
                description, isMultitoken);
          }
      }
  }

  /**
   * Add an advanced parameter (can be a multitoken or not).
   *
   * @param applicationParameter: advanced parameter.
   * @param isMultitoken: does the parameter contain multiple values?
   * @param configOption: is it a configuration parameter (inside the file)?
   * @param hideOption: is option displaying (when we call the application with
   *        "help" parameter)?
   * @param genericOption: is it a generic Option (like --help, --version...)?
   */
  template<typename T, typename D>
  void addParameter(ApplicationAdvancedParameter<T, D>& applicationParameter,
                    bool isMultitoken,
                    bool configOption = true,
                    bool hideOption = true,
                    bool genericOption = false)
  {
    addAdvancedParameter(applicationParameter.name(),
                         applicationParameter.valuePtr(),
                         applicationParameter.description(),
                         isMultitoken, configOption, hideOption, genericOption);
  }


private:
  template<typename T>
  void addParameterTo(boost::program_options::options_description& od,
                      const std::string& name, T* value, const T& defaultValue,
                      const std::string& description)
  {
    od.add_options()(name.c_str(),
            boost::program_options::value<T>(value)->default_value(defaultValue),
            description.c_str());
  }

  template<typename T>
  void addParameterTo(boost::program_options::options_description& od,
                      const std::string& name, T* value,
                      const std::string& description, bool isMultitoken)
  {
    if (!isMultitoken)
      {
        od.add_options()(name.c_str(), boost::program_options::value(value), description.c_str());
      }
    else
      {
        od.add_options()(name.c_str(), boost::program_options::value(value)->multitoken(), description.c_str());
      }
  }
};

template<class T, typename U>
void validate(boost::any & v, const std::vector<std::string> & values, T* paramPtr, U type)
{
  T tmpPtr;
  tmpPtr.validate(v, values);
}

template<typename U>
void validate(boost::any & v, const std::vector<std::string> & values, Double* paramPtr, U type)
{
  //does nothing.
}

#endif /* BCPARAMETERMANAGER_HPP_ */
