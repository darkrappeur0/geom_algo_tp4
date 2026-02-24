#pragma once
//lissage hpp
#include "types.hpp"

namespace geomAlgoLib
{
    Mesh lissage_1_iter(const Mesh & myMesh);
    Mesh lissage_n_iter(const Mesh & myMesh, int n);
}