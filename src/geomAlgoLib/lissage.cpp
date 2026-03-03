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

    Mesh lissage_de_taubin(const Mesh& myMesh, double lambda, double mu)
{
    Mesh temp_mesh = myMesh;

    std::vector<Kernel::Point_3> positions_lambda;
    std::vector<Kernel::Point_3> positions_mu;

    // ---------------------
    // 1ère passe (lambda)
    // ---------------------
    for(auto v = temp_mesh.vertices_begin(); v != temp_mesh.vertices_end(); ++v)
    {
        CGAL::Vector_3<Kernel> laplacian(0,0,0);
        int degree = v->vertex_degree();

        auto h = v->vertex_begin();
        for(int i = 0; i < degree; ++i, ++h)
        {
            auto u = h->opposite()->vertex();
            laplacian += (u->point() - v->point());
        }

        laplacian /= degree;

        positions_lambda.push_back(v->point() + lambda * laplacian);
    }

    auto v_it = temp_mesh.vertices_begin();
    for(size_t i = 0; i < positions_lambda.size(); ++i, ++v_it)
        v_it->point() = positions_lambda[i];


    for(auto v = temp_mesh.vertices_begin(); v != temp_mesh.vertices_end(); ++v)
    {
        CGAL::Vector_3<Kernel> laplacian(0,0,0);
        int degree = v->vertex_degree();

        auto h = v->vertex_begin();
        for(int i = 0; i < degree; ++i, ++h)
        {
            auto u = h->opposite()->vertex();
            laplacian += (u->point() - v->point());
        }

        laplacian /= degree;

        positions_mu.push_back(v->point() + mu * laplacian);
    }

    // appliquer mu
    v_it = temp_mesh.vertices_begin();
    for(size_t i = 0; i < positions_mu.size(); ++i, ++v_it)
        v_it->point() = positions_mu[i];

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
                double weight = 1.0 / (std::sqrt(vector.squared_length()));
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

    
    Mesh lissage_1_iter_ponderer_contangentes(const Mesh & myMesh){
        Mesh temp_mesh = myMesh;
        std::vector<Kernel::Point_3> new_positions;
        for(auto v = myMesh.vertices_begin(); v != myMesh.vertices_end(); ++v){
            auto c = v->vertex_degree();
            auto h = v->vertex_begin();
            CGAL::Vector_3<Kernel> vector_total (0.0,0.0,0.0); 
            double weight_sum = 0.0;
            while (c!=0){
                auto h2 = h->opposite();
                auto u = h2->vertex();
                
                auto hnext = h->next();
                auto p = hnext->vertex(); 

                auto vectorv1aij =  v->point() - p->point();
                auto vectorv2aij =  u->point() - p->point();
                
                double norme_v1 = std::sqrt(vectorv1aij.squared_length());
                double norme_v2 = std::sqrt(vectorv2aij.squared_length());
                double cotanaij = CGAL::scalar_product(vectorv1aij, vectorv2aij) /std::sqrt(CGAL::cross_product(vectorv1aij, vectorv2aij).squared_length());
                
                
                
                auto h2next = h2->next();
                auto p2 = h2next->vertex();
                auto vectorv1bij =  v->point() - p2->point();
                auto vectorv2bij =  u->point() - p2->point();
                
                double cotanbij = CGAL::scalar_product(vectorv1bij, vectorv2bij) /std::sqrt(CGAL::cross_product(vectorv1bij, vectorv2bij).squared_length());
                double weight = 0.5 * (cotanaij + cotanbij);

                
                auto vector = u->point() - v->point();
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

    Mesh lissage_n_iter_ponderer_contangentes(const Mesh & myMesh, int n){
        Mesh temp_mesh = myMesh;
        for(int i = 0;i<n;i++){
            temp_mesh = lissage_1_iter_ponderer_contangentes(temp_mesh);
        }
        return temp_mesh;

    }
        




}

