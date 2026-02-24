#include <io.hpp>
#include <example.hpp>
#include <lissage.hpp>

#include <iostream>
#include <string>

int main(int argc, char *argv[]){

    std::cout << "Hello !" << std::endl;

    if(argc < 2){
        throw std::invalid_argument("This program expects at least 1 argument (path to a mesh).");
    }

    const std::string meshPath = std::string{argv[1]};
    
    geomAlgoLib::Mesh myMesh;

    geomAlgoLib::readOFF(meshPath, myMesh);
    geomAlgoLib::remove_degenerate_faces(myMesh);

    auto genus = geomAlgoLib::computeGenus(myMesh);
    std::cout << "The Genus of [" << meshPath << "] is = " << std::to_string(genus) << std::endl;

    geomAlgoLib::writeOFF(myMesh,"output.off");

    /*
    geomAlgoLib::FacetDoubleMap tab = geomAlgoLib::aire_calcul(myMesh);
	std::cout << "apres le calcul de l'aire au carrée" << std::endl;
    unsigned c = 0;
    for ( auto i = tab.begin();i!=tab.end();++i){
        c=c+1;
        std::cout << "aire a la " << c << "-eme valeur : " <<i->second << std::endl;
        CGAL::Vector_3<geomAlgoLib::Kernel> v = geomAlgoLib::calcul_angle(i->first);
        std::cout << "voici la valeur de v " << v << std::endl;
        
        
    }
    geomAlgoLib::StoreinFiles(myMesh,tab,"test.off");

    geomAlgoLib::FacetDoubleMap angle = geomAlgoLib::angle_min_calcul(myMesh);

    c=0;

    for(auto i = angle.begin();i!=angle.end();++i){
        c=c+1;
        std::cout << "plus petit angle en degrés de la " << c << "-eme face : " <<i->second << " et son aire : " << tab[i->first] << std::endl;
    }

    geomAlgoLib::FacetStringMap res = geomAlgoLib::etiquettage(tab);
    
    c = 0;
    for(auto i = res.begin();i!=res.end();++i){
        c=c+1;
        std::cout << "etiquette de la " << c << "-eme face : " <<i->second << " et son aire : " << tab[i->first] << std::endl;
    }
    

    
    
    geomAlgoLib::StoreinFiles_labels(myMesh,res,tab,"etiquettes.off");
    
    //std::cout << "tableau des aires" << tab << std::endl;
    */


    auto newMesh = geomAlgoLib::lissage_1_iter(myMesh);

    geomAlgoLib::writeOFF(newMesh,"output2.off");


    auto newMesh2 = geomAlgoLib::lissage_n_iter(myMesh,4);

    geomAlgoLib::writeOFF(newMesh2,"output3.off");






    std::cout << "The end..." << std::endl;
    return 0;
}