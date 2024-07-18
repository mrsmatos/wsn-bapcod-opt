#ifndef VRPSTW_LOADER_H
#define VRPSTW_LOADER_H

#include <string>

namespace vrpstw
{
	class Data;
	class Parameters;

	class Loader
	{
	public:
		Loader();

		bool loadData(const std::string & file_name);
		bool loadParameters(const std::string & file_name, int argc, char* argv[]);

	private:
		Data & data;
		Parameters & parameters;

		bool loadVRPTWFile(std::ifstream & ifs);
	};
}

#endif
