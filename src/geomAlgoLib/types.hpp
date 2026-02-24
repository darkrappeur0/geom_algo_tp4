#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Kernel/global_functions_3.h>

#include <map>

namespace geomAlgoLib{

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Polyhedron_3<Kernel>;

using FacetCstIt = Mesh::Facet_const_iterator;
using VertexCstIt = Mesh::Vertex_const_iterator;
using HalfedgeCstIt = Mesh::Halfedge_const_iterator;
using HalfedgeFacetCstCirc = Mesh::Halfedge_around_facet_const_circulator;

using FacetDoubleMap = std::map<Mesh::Facet_const_handle, double>;
using FacetIntMap = std::map<Mesh::Facet_const_handle, int>;
using FacetStringMap = std::map<Mesh::Facet_const_handle, std::string>;

}