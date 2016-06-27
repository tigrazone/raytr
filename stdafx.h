//
// stdafx.h: Cache for common headers
// Created by dc on 5/7/16.
//

#ifndef RAYTR_STDAFX_H
#define RAYTR_STDAFX_H

// STL

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <atomic>
#include <set>
#include <thread>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <parallel/algorithm>
#include <parallel/numeric>

namespace GP = __gnu_parallel;

using std::tie;
using std::pair;
using std::tuple;
using std::make_tuple;
using std::make_pair;
using std::vector;
using std::string;
using std::cerr;
using std::clog;
using std::min;
using std::max;
using std::swap;
using std::endl;
using std::shared_ptr;
using std::unique_ptr;
using std::make_shared;

// JSON

#include <json.hpp>
using nlohmann::json;

#define J_AT_F(J, k) J.at(k).get<double>()
#define J_AT_I(J, k) J.at(k).get<int>()
#define J_AT_F_NZ(j, k, def) (j.count(k) ? j.at(k).get<double>(): def)
#define J_AT_I_NZ(j, k, def) (j.count(k) ? j.at(k).get<int>(): def)

#include <Eigen/Eigen>

#endif //RAYTR_STDAFX_H
