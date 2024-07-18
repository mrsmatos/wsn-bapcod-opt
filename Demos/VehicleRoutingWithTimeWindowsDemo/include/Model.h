#ifndef VRPSTW_MODEL_H
#define VRPSTW_MODEL_H

#include <bcModelPointerC.hpp>
#include "InputUser.h"

namespace vrpstw
{
	class Model : public BcModel, InputUser
	{
	public:
		Model(const BcInitialisation& bc_init);
		virtual ~Model() {}
	};
}

#endif
