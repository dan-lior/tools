#pragma once

#include "common_defs.h"
#include "voxel_box.h"

std::vector<GridWithAttribute<GridIso<3>, Vector<3>>> read_nondicom(const std::string& nondicom_directory);
