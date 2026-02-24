#pragma once

#include "types.hpp"

namespace geomAlgoLib
{
     /*
        Toutes les fonctions interne apparaisse ici si vous voulez les tester dans le main
        indépendament de leur super fonctions si vous voulez les tester.

        Pour celle qui sont commenter oubliez pas de décommenter le code associé dans measures.cpp
     
     
     */
    
    //double calcul_air_triangle(const FacetCstIt &f);
    //double calcul_quad(const FacetCstIt &i);
    double calcul_air_face(const FacetCstIt &f);
    FacetDoubleMap aire_calcul(const Mesh &mesh);


    CGAL::Vector_3<Kernel> calcul_normal(const FacetCstIt &f);
    CGAL::Vector_3<Kernel> calcul_angle(const FacetCstIt &f);
    
    
    bool est_zone_marche(const FacetCstIt& f);
    bool est_obstacle(const FacetCstIt& f);
    bool est_a_portee(const FacetCstIt& f);
    FacetStringMap etiquettage(const FacetDoubleMap &f);


    double calcul_angle_min_face(const FacetCstIt& f);
    FacetDoubleMap angle_min_calcul(const Mesh &mesh);
    void remove_degenerate_faces(Mesh &mesh);

}