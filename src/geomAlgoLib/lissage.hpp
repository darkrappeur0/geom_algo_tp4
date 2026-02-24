#pragma once
//lissage hpp
#include "types.hpp"

namespace geomAlgoLib
{
    Mesh lissage_1_iter(const Mesh & myMesh);
    Mesh lissage_n_iter(const Mesh & myMesh, int n);

    Mesh facteur_de_diffusion(const Mesh & myMesh,double lambda);
    Mesh facteur_de_diffusion_n_iter(const Mesh & myMesh, int n,double lambda);

    Mesh lissage_de_taubin(const Mesh & myMesh,double lambda,double mu);
    Mesh lissage_de_taubin_n_iter(const Mesh & myMesh, int n,double lambda,double mu);

    Mesh lissage_1_iter_ponderer(const Mesh & myMesh);
    Mesh lissage_n_iter_ponderer(const Mesh & myMesh, int n);
}