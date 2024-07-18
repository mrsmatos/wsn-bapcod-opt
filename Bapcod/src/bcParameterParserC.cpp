/**
 *
 * This file bcParameterParserC.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcParameterParserC.hpp"
#include "bcApplicationException.hpp"
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>

using namespace std;
namespace po = boost::program_options;

ParameterParser::ParameterParser(
    std::string parameterFileName) : _parameterFileName(parameterFileName), _isAlreadyParsed(false)
{
}

ParameterParser::~ParameterParser()
{
}

void ParameterParser::parse(int argc, char** argv)
{
  if (!_isAlreadyParsed)
    {
      //Application parameters
      po::options_description generic("Generic options");

      generic.add_options()("usage,u", "Produce help message")
                           ("applicationParameters,a",
                            po::value<string>(&_parameterFileName)->default_value(_parameterFileName),
                            "Path to the application specific parameters file.");

      _cmdLineOptions.add(generic);
      _visibleOptions.add(generic);

      po::variables_map vm;

      if (argc > 0) {
          try {
              po::store(po::command_line_parser(argc, argv).options(_cmdLineOptions).allow_unregistered().run(), vm);
              po::notify(vm);


#if(BOOST_VERSION > 104000)
          } catch (const boost::program_options::multiple_occurrences &e) {
              throw ApplicationException("Can't parse the file: " + std::string(e.what()) + ". Option: "
                                         + std::string(e.get_option_name()));
#endif
          } catch (const exception &e) {
              throw ApplicationException("Can't parse the file: " + std::string(e.what()));
          }
      }

      if(vm.count("usage"))
        {
          cout << _visibleOptions << endl;
          exit(0);
        }

      ifstream ifs(_parameterFileName.c_str());
      if (!ifs)
        {
          throw ApplicationException("Can not open application config file: " + _parameterFileName,
                                     false);
        }
      else
        {
          try{
              po::store(po::parse_config_file(ifs, _configFileOptions, false), vm);
              po::notify(vm);
#if(BOOST_VERSION > 104000)
          } catch (const boost::program_options::multiple_occurrences& e) {
              throw ApplicationException("Can't parse the file: " + std::string(e.what()) + ". Option: "
                                         + std::string(e.get_option_name()));
#endif
          }catch (const exception& e) {
              throw ApplicationException("Can't parse the file: " + std::string(e.what()));
          }
        }
      _isAlreadyParsed = true;
    }
  else {
      throw ApplicationException("The parameter file is already parsed: " + _parameterFileName);
  }
}


void ParameterParser::parseFile(std::string parameterFileName){

    boost::program_options::options_description newCmdLineOptions;
    po::options_description generic("Generic options");

    generic.add_options()("usage,u", "Produce help message")
                         ("applicationParameters,a",
                          po::value<string>(&parameterFileName)->default_value(parameterFileName),
                          "Path to the application specific parameters file.");

    newCmdLineOptions.add(generic);
    ifstream ifs(parameterFileName.c_str());

    po::variables_map vm;


    try{
        po::store(po::command_line_parser(0,NULL).options(newCmdLineOptions).allow_unregistered().run(),
                  vm);
        po::notify(vm);


#if(BOOST_VERSION > 104000)
    } catch (const boost::program_options::multiple_occurrences& e) {
        throw ApplicationException("Can't parse the file: " + std::string(e.what()) + ". Option: "
                                   + std::string(e.get_option_name()));
#endif
    }catch (const exception& e) {
        throw ApplicationException("Can't parse the file: " + std::string(e.what()));
    }

    try{
        po::store(po::parse_config_file(ifs, _configFileOptions, false), vm);

        po::notify(vm);
#if(BOOST_VERSION > 104000)
    } catch (const boost::program_options::multiple_occurrences& e) {
        throw ApplicationException("Can't parse the file: " + std::string(e.what()) + ". Option: "
                                    + std::string(e.get_option_name()));
#endif
    }catch (const exception& e) {
        throw ApplicationException("Can't parse the file: " + std::string(e.what()));
    }
}


void ParameterParser::setParameterFileName(const std::string & fileName)
{
    _parameterFileName = fileName;
}

std::string ParameterParser::getParameterFileName() const
{
  return _parameterFileName;
}
