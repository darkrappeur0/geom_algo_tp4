#include "lissage.hpp"
//lissage cpp
namespace geomAlgoLib
{

    Mesh lissage_1_iter(const Mesh & myMesh){
        Mesh temp_mesh = myMesh;
        std::vector<Kernel::Point_3> new_positions;
        for(auto v = myMesh.vertices_begin(); v != myMesh.vertices_end(); ++v){
            auto c = v->vertex_degree();
            auto h = v->vertex_begin();
            CGAL::Vector_3<Kernel> vector_total (0.0,0.0,0.0); 
            while (c!=0){
                auto u = h->opposite()->vertex();
                auto vector = u->point() - v->point();
                vector_total = vector_total + vector;   
                ++h;
                c=c-1;
            }
            vector_total = vector_total/v->vertex_degree();
            double lambda = 1.0;

             Kernel::Point_3 new_p =v->point() + lambda * vector_total;

            new_positions.push_back(new_p);
        }
        auto v2 = temp_mesh.vertices_begin();
        for(size_t i = 0; i < new_positions.size(); ++i, ++v2)
        {
            v2->point() = new_positions[i];
        }
        return temp_mesh;
    } 

    Mesh lissage_n_iter(const Mesh & myMesh, int n){
        Mesh temp_mesh = myMesh;
        for(int i = 0;i<n;i++){
            temp_mesh = lissage_1_iter(temp_mesh);
        }
        return temp_mesh;

    }

    Mesh facteur_de_diffusion(const Mesh & myMesh,double lambda){
        Mesh temp_mesh = myMesh;
        std::vector<Kernel::Point_3> new_positions;
        for(auto v = myMesh.vertices_begin(); v != myMesh.vertices_end(); ++v){
            auto c = v->vertex_degree();
            auto h = v->vertex_begin();
            CGAL::Vector_3<Kernel> vector_total (0.0,0.0,0.0); 
            while (c!=0){
                auto u = h->opposite()->vertex();
                auto vector = u->point() - v->point();
                vector_total = vector_total + vector;   // faire le calcul de la position du nouveau point.
                ++h;
                c=c-1;
            }
            vector_total = vector_total/v->vertex_degree();
            

            Kernel::Point_3 new_p =v->point() + lambda * vector_total;

            new_positions.push_back(new_p);
        }
        auto v2 = temp_mesh.vertices_begin();
        for(size_t i = 0; i < new_positions.size(); ++i, ++v2)
        {
            v2->point() = new_positions[i];
        }
        return temp_mesh;
    } 

    Mesh facteur_de_diffusion_n_iter(const Mesh & myMesh, int n,double lambda){
        Mesh temp_mesh = myMesh;
        for(int i = 0;i<n;i++){
            temp_mesh = facteur_de_diffusion(temp_mesh,lambda);
        }
        return temp_mesh;

    }

    Mesh lissage_de_taubin(const Mesh & myMesh,double lambda,double mu){
        Mesh temp_mesh = myMesh;
        std::vector<Kernel::Point_3> new_positions;
        for(auto v = myMesh.vertices_begin(); v != myMesh.vertices_end(); ++v){
            auto c = v->vertex_degree();
            auto h = v->vertex_begin();
            CGAL::Vector_3<Kernel> vector_total (0.0,0.0,0.0); 
            while (c!=0){
                auto u = h->opposite()->vertex();
                auto vector = u->point() - v->point();
                vector_total = vector_total + vector;   // faire le calcul de la position du nouveau point.
                ++h;
                c=c-1;
            }
            vector_total = vector_total/v->vertex_degree();
            

            Kernel::Point_3 new_p =v->point() + lambda * vector_total;
            new_p =new_p + mu * vector_total;

            new_positions.push_back(new_p);
        }
        auto v2 = temp_mesh.vertices_begin();
        for(size_t i = 0; i < new_positions.size(); ++i, ++v2)
        {
            v2->point() = new_positions[i];
        }
        return temp_mesh;
    } 

    Mesh lissage_de_taubin_n_iter(const Mesh & myMesh, int n,double lambda,double mu){
        Mesh temp_mesh = myMesh;
        for(int i = 0;i<n;i++){
            temp_mesh = lissage_de_taubin(temp_mesh,lambda,mu);
        }
        return temp_mesh;

    }



    Mesh lissage_1_iter_ponderer(const Mesh & myMesh){
        Mesh temp_mesh = myMesh;
        std::vector<Kernel::Point_3> new_positions;
        for(auto v = myMesh.vertices_begin(); v != myMesh.vertices_end(); ++v){
            auto c = v->vertex_degree();
            auto h = v->vertex_begin();
            CGAL::Vector_3<Kernel> vector_total (0.0,0.0,0.0); 
            double weight_sum = 0.0;
            while (c!=0){
                auto u = h->opposite()->vertex();
                auto vector = u->point() - v->point();
                double weight = 1.0 / (std::sqrt(vector.squared_length()) + 1e-8);
                vector_total = vector_total + vector*weight;   
                weight_sum += weight;
                ++h;
                c=c-1;
            }
            if(weight_sum > 0) vector_total /= weight_sum;
            

            Kernel::Point_3 new_p =v->point() + vector_total;

            new_positions.push_back(new_p);
        }
        auto v2 = temp_mesh.vertices_begin();
        for(size_t i = 0; i < new_positions.size(); ++i, ++v2)
        {
            v2->point() = new_positions[i];
        }
        return temp_mesh;
    } 

    Mesh lissage_n_iter_ponderer(const Mesh & myMesh, int n){
        Mesh temp_mesh = myMesh;
        for(int i = 0;i<n;i++){
            temp_mesh = lissage_1_iter_ponderer(temp_mesh);
        }
        return temp_mesh;

    }


}

