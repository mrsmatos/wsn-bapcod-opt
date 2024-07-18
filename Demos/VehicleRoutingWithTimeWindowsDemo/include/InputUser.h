#ifndef VRPSTW_INPUTUSER_H
#define VRPSTW_INPUTUSER_H

#include "Data.h"
#include "Parameters.h"

namespace vrpstw
{
	class InputUser
	{
	protected:
		InputUser() : data(Data::getInstance()), params(Parameters::getInstance()) {}

		const Data & data;
		const Parameters & params;
	};
}

#endif
