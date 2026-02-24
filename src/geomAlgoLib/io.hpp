#pragma once

#include "types.hpp"

#include <iostream>
#include <fstream>

namespace geomAlgoLib
{

/// Read an OFF file and store the mesh in "mesh".
/// Returns true if it was loaded successfully.
bool readOFF(const std::string& filePath, Mesh& mesh);

/// Write a mesh at location "filePath"
void writeOFF(const Mesh& mesh, const std::string& filePath);

void StoreinFiles(const Mesh& mesh,const FacetDoubleMap & calcul, const std::string& filePath );

double calcul_couleur_mesh(const double& second);

void StoreinFiles_labels(const Mesh& mesh,const FacetStringMap & tab1,const FacetDoubleMap & tab2 ,const std::string& filePath );

}