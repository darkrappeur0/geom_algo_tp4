#include "measures.hpp"
#include <cmath>
#include <exception>
namespace geomAlgoLib
{
    /*
    
    
    Code non utiliser, présent si vous voulez voir le précédent code utilisé
    A décommenter si vous voulez les tester
    
    */
    /*
    double calcul_air_triangle(const FacetCstIt &f) {
        HalfedgeCstIt edges_1 = f->halfedge();
        HalfedgeCstIt edges_2 = edges_1->next();
        HalfedgeCstIt edges_3 = edges_2->next();
        
        double r =0.0;
        const auto& p = edges_1->vertex()->point();
        const auto& q = edges_2->vertex()->point();
        const auto& s = edges_3->vertex()->point();
        CGAL::Vector_3<Kernel> pq = q - p;
        CGAL::Vector_3<Kernel> ps = s - p;
        std::cout << "Je suis un triangle" << std::endl;
        
        if(!CGAL::collinear(p,q,s) ){
            r =std::sqrt(CGAL::cross_product(pq, ps).squared_length()) * 0.5;
        }
        return r; 
    }
        */

    /*    
    double calcul_quad(const FacetCstIt &i) {
        HalfedgeCstIt edges_1 = i->halfedge();
        HalfedgeCstIt edges_2 = edges_1->next();
        HalfedgeCstIt edges_3 = edges_2->next();
        HalfedgeCstIt edges_4 = edges_3->next();

        double r1 = 0.0; 
        double r2 = 0.0 ;
        std::cout << "Je suis un carré" << std::endl;
        try{
            if (!CGAL::collinear(edges_1->vertex()->point(),edges_2->vertex()->point(),edges_3->vertex()->point())){
                r1 = CGAL::squared_area(
                edges_1->vertex()->point(),
                edges_2->vertex()->point(),
                edges_3->vertex()->point()
                );
            }
        
            if (!CGAL::collinear(edges_3->vertex()->point(),edges_4->vertex()->point(),edges_1->vertex()->point())){
                r2 = CGAL::squared_area(
                edges_3->vertex()->point(),
                edges_4->vertex()->point(),
                edges_1->vertex()->point()
                );
            }
        }
        catch (const CGAL::Precondition_exception & e){
            r1=0.0;
            r2=0.0;
        }
        r1 = std::sqrt(r1);
        r2 = std::sqrt(r2);
        return r1+r2;
    }
        */
    
    /*
    FacetDoubleMap aire_calcul(const Mesh &mesh){
        
        FacetDoubleMap tab;
        for (FacetCstIt i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	    {

            if (i->is_triangle()){
                
                tab[i] = calcul_air_triangle(i);
                
                
            }
            else{
                if (i->is_quad()){
                    tab[i] = calcul_quad(i);
                }
                else{
                    throw std::runtime_error("Face ni triangle ni quadrilatère");
                }
            }
           
		        
	    }
        return tab;
   } 
        */


    double calcul_air_face(const FacetCstIt &f)
{
    HalfedgeCstIt h = f->halfedge();
    HalfedgeCstIt start = h;

    // Point pivot
    const auto& p0 = h->vertex()->point();

    h = h->next();

    double aire_totale = 0.0;

    while(h->next() != start)
    {
        const auto& p1 = h->vertex()->point();
        const auto& p2 = h->next()->vertex()->point();

        if(!CGAL::collinear(p0, p1, p2))
        {
            CGAL::Vector_3<Kernel> v1 = p1 - p0;
            CGAL::Vector_3<Kernel> v2 = p2 - p0;

            double aire_triangle =
                std::sqrt(CGAL::cross_product(v1, v2).squared_length()) * 0.5;

            aire_totale += aire_triangle;
        }

        h = h->next();
    }

    return aire_totale;
}


FacetDoubleMap aire_calcul(const Mesh &mesh){
        
        FacetDoubleMap tab;
        for (FacetCstIt i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	    {

            tab[i] = calcul_air_face(i);
		        
	    }
        return tab;
   } 


   CGAL::Vector_3<Kernel> calcul_normal(const FacetCstIt &f){
        HalfedgeCstIt edges_1 = f->halfedge();
        HalfedgeCstIt edges_2 = edges_1->next();
        HalfedgeCstIt edges_3 = edges_2->next();
        return CGAL::normal(edges_1->vertex()->point(),edges_2->vertex()->point(),edges_3->vertex()->point());
   }


   CGAL::Vector_3<Kernel> calcul_angle(const FacetCstIt &f){
    CGAL::Vector_3<Kernel> n = calcul_normal(f);
    n = n / std::sqrt(n.squared_length());
    double theta_x = std::acos(n.x());
    double theta_y = std::acos(n.y());
    double theta_z = std::acos(n.z());

    double rad_to_deg = 180.0 / 3.14;

    theta_x *= rad_to_deg;
    theta_y *= rad_to_deg;
    theta_z *= rad_to_deg;
    CGAL::Vector_3<Kernel> v{Kernel::FT(theta_x), Kernel::FT(theta_y), Kernel::FT(theta_z)};

    return v;
   }



bool est_zone_marche(const FacetCstIt& f)
{
    CGAL::Vector_3<Kernel> angles = calcul_angle(f);
    
    
    if (angles.x() == 0 && angles.y() == 0 && angles.z() == 0)
        return false;
    
    double theta_z_deg = CGAL::to_double(angles.z());
    
    
    return (theta_z_deg <= 25.84);
}

  bool est_obstacle(const FacetCstIt& f)
{
    CGAL::Vector_3<Kernel> angles = calcul_angle(f);
    
    
    if (angles.x() == 0 && angles.y() == 0 && angles.z() == 0)
        return false;
    
    double theta_z_deg = CGAL::to_double(angles.z());
    
    
    return (theta_z_deg > 25.84 );
}


    bool est_a_portee(const FacetCstIt& f)
    {
    HalfedgeCstIt h = f->halfedge();

    Kernel::Point_3 p = h->vertex()->point();
    Kernel::Point_3 q = h->next()->vertex()->point();
    Kernel::Point_3 r = h->next()->next()->vertex()->point();

    Kernel::Point_3 bary(
        (p.x() + q.x() + r.x())/3.0,
        (p.y() + q.y() + r.y())/3.0,
        (p.z() + q.z() + r.z())/3.0
    );

    Kernel::Point_3 avatar(3.0, 3.0, 0.0);

    double dist2 = CGAL::squared_distance(bary, avatar);

    return dist2 < 4.0; // distance < 2
    }


   FacetStringMap etiquettage(const FacetDoubleMap &tab){
    FacetStringMap res;
    for(auto i = tab.begin(); i != tab.end(); ++i){
        if (i->second > 0.0005){
            res[i->first] = "Grande Face";
        }
        else{
            res[i->first] = "Petite Face";
        }
        
        CGAL::Vector_3<Kernel> n = calcul_normal(i->first);
        n = n / std::sqrt(n.squared_length());
        if( est_zone_marche(i->first) & (res[i->first] == "Grande Face")){
            res[i->first] = "orientées vers le haut";
        }
        if( est_obstacle(i->first) ){
            res[i->first] = "obstacle";
        }
        if( est_a_portee(i->first) ){
            res[i->first] = "a portee";
        }
        
    }
    return res;
   }


   double calcul_angle_min_face(const FacetCstIt& f)
{
    HalfedgeCstIt h = f->halfedge();
    HalfedgeCstIt start = h;

    double angle_min = std::numeric_limits<double>::max();

    do
    {
        // Points
        const auto& p_prev = h->prev()->vertex()->point();
        const auto& p      = h->vertex()->point();
        const auto& p_next = h->next()->vertex()->point();

        // Vecteurs des deux arêtes adjacentes au sommet p
        CGAL::Vector_3<Kernel> v1 = p_prev - p;
        CGAL::Vector_3<Kernel> v2 = p_next - p;

        double norm1 = std::sqrt(v1.squared_length());
        double norm2 = std::sqrt(v2.squared_length());

        if(norm1 > 0.0 && norm2 > 0.0)
        {
            double cos_theta = (v1 * v2) / (norm1 * norm2);

            
            cos_theta = std::max(-1.0, std::min(1.0, cos_theta));

            double angle = std::acos(cos_theta);

            angle_min = std::min(angle_min, angle);
        }

        h = h->next();

    } while(h != start);

    double rad_to_deg = 180.0 / CGAL_PI;
    return angle_min * rad_to_deg;

}


FacetDoubleMap angle_min_calcul(const Mesh &mesh)
{
    FacetDoubleMap tab;

    for (FacetCstIt i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
    {
        tab[i] = calcul_angle_min_face(i);
    }

    return tab;
}



void remove_degenerate_faces(Mesh &mesh)
{

    std::vector<Mesh::Facet_iterator> to_remove;
    
    for (Mesh::Facet_iterator f = mesh.facets_begin(); f != mesh.facets_end(); ++f)
    {
        Mesh::Halfedge_handle h = f->halfedge();
        
        if (h->next()->next()->next() == h)
        {
            const auto& p0 = h->vertex()->point();
            const auto& p1 = h->next()->vertex()->point();
            const auto& p2 = h->next()->next()->vertex()->point();
            
            if (CGAL::collinear(p0, p1, p2))
            {
                to_remove.push_back(f);
            }
        }
    }
    
    std::cout << "Suppression de " << to_remove.size() << " faces dégénérées..." << std::endl;
    
    for (auto f : to_remove)
    {
        mesh.erase_facet(f->halfedge());
    }
}


}