#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <ios>
#include <streambuf>
#include <map>
#include <set>
#include <array>
#include <random>
#include <type_traits>
#include <exception>
#include <algorithm>
#include <chrono>
#include <thread>
#include <atomic>
#include <mutex>
#include <tuple>
#include <omp.h>
#include <time.h>
#include <cmath>
#include <cassert>

#define EIGEN_DEFAULT_TO_COLUMN_MAJOR
#include "./include/Eigen_3.3.1/Eigen/Core"
#include "./include/Eigen_3.3.1/Eigen/Dense"

#include "./include/gsl-2.5/gsl/gsl_vector.h"
#include "./include/gsl-2.5/gsl/gsl_multimin.h"

#include "EigenDummyTensor.h"
#include "EigenMultiDimArray.h"
