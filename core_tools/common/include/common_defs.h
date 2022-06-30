#pragma once

#include <math.h>
#include <cassert>
#include <string>
#include <array>
#include <vector>
#include <stack>
#include <map>
#include <set>
#include <cstddef> // for uint64_t
#include <cmath>
#include <cstdint> //?

#include <iostream>
#include <fstream>
//#include <chrono>
#include <algorithm>
#include <limits>
#include <cassert>
#include <random>
#include <filesystem>
//#include <optional>
#include <Eigen/Dense>


template<uint64_t dim>
using Matrix = Eigen::Matrix<double, dim, dim>;

template<uint64_t dim_rows, uint64_t dim_cols>
using MatrixRect = Eigen::Matrix<double, dim_rows, dim_cols>;

template<uint64_t dim>
using Vector = Eigen::Matrix<double, dim, 1>;

template<uint64_t dim>
using Index = Eigen::Matrix<int64_t, dim, 1>;

// cast each component to double
template<uint64_t dim>
Vector<dim> index_to_vector(const Index<dim>& index);

// performs a floor operation on each component
template<uint64_t dim>
Index<dim> vector_to_index(const Vector<dim>& vec);

#include "common_defs.hpp"