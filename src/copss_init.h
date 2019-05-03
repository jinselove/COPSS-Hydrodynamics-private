// Copyright (C) 2015-2016 Jiyuan Li, Xikai Jiang

// This code is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.


// This code is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.


// You should have received a copy of the GNU General Public
// License along with this code; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#pragma once

// libMesh
#include "libmesh/libmesh.h"

/**
 * Initialization object for any Copssd application
 *
 * This object must be created in the main() of any Copss application so
 * everything is properly initialized and finalized.
 */
namespace libMesh {
class CopssInit : public LibMeshInit {
public:

  CopssInit(int    argc,
            char **argv);

  virtual ~CopssInit() = default;
};
}
