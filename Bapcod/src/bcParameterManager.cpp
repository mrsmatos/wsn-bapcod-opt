/**
 *
 * This file bcParameterManager.cpp is a part of BaPCod - a generic Branch-And-Price Code.
 * https://bapcod.math.u-bordeaux.fr
 *
 * Inria Bordeaux, All Rights Reserved. 
 * See License.pdf file for license terms.
 * 
 * This file is available only for academic use. 
 * For other uses, please contact stip-bso@inria.fr.
 *
 **/

#include "bcParameterManager.hpp"
#include "bcGlobalException.hpp"

#include <boost/program_options.hpp>
#include "bcBapcodVersion.hpp"

#include <stdexcept>

using namespace std;
namespace po = boost::program_options;

ParameterManager::ParameterManager() :
    _configFileOptions("Config file options"),
    _cmdLineOptions("Command line options"),
    _visibleOptions("Visible option for the help command"),
    _genericOptions("Help command") {
}

ParameterManager::~ParameterManager() {
}

void ParameterManager::parse(int argc, char** argv, std::string& filename)
{
    _genericOptions.add_options()
            ("version,v", "Print version string")
            ("help,h", "Produce help message");
    _cmdLineOptions.add(_genericOptions);
    _visibleOptions.add(_genericOptions);

    po::variables_map vm;

    std::string filenameCopy = filename; //added by Issam.

    if (argc > 0) {
        try {
            po::store(po::command_line_parser(argc, argv).options(_cmdLineOptions).allow_unregistered().run(),
                      vm);
            po::notify(vm);
            if (filename == "NOT_SPECIFIED") {
                filename.assign(filenameCopy);
            }
        } catch (const boost::program_options::multiple_occurrences &e) {
            throw GlobalException("Can't parse the file: " + filename + " because: " +
                                  std::string(e.what()) + ". Option: " + std::string(e.get_option_name()),
                                  true);
        }
    }

    if (vm.count("help")) {
        cout << _visibleOptions << endl;
        exit(EXIT_SUCCESS);
    }

    if (vm.count("version")) {
        cout << "Bapcod Version: " << BAPCOD_VERSION << endl;
        exit(EXIT_SUCCESS);
    }

    ifstream ifs(filename.c_str());
    //If the file is found (and opened)

    if (!ifs) {
        throw GlobalException("Can not open bapcod config file: " + filename, true);
    }
    
    if (ifs) {
        //We parsed it.
        try {
            //Issam: unrecognized options should not be ignored.
//            po::store(po::parse_config_file(ifs, _configFileOptions, true), vm);

            po::store(po::parse_config_file(ifs, _configFileOptions, false), vm);
          
            po::notify(vm);
            //if block added by Issam.
            if (filename == "NOT_SPECIFIED")
            {
                filename.assign(filenameCopy);
            }
        } catch (const boost::program_options::multiple_occurrences& e) {
            throw GlobalException("Can't parse the file: " + filename + " because: " +
                                  std::string(e.what()) + ". Option: " + std::string(e.get_option_name()),
                                  true);
        } catch (const boost::program_options::error & ee) {
            throw GlobalException("Can't parse the file: " + filename + " because: " + std::string(ee.what()),
                                  true);
        }
    }
    else {
        cout << "Config file not found " << filename << endl;
    }
}
