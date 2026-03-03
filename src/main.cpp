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

    geomAlgoLib::writeOFF(myMesh,"meshdebase.off");

    /*
    
        lissage laplacien 1 iter

    */    
    auto newMesh = geomAlgoLib::lissage_1_iter(myMesh);

    geomAlgoLib::writeOFF(newMesh,"output2.off");

    /*
    
        lissage laplacien n iter
    
    */
    auto newMesh2 = geomAlgoLib::lissage_n_iter(myMesh,4);

    geomAlgoLib::writeOFF(newMesh2,"output3.off");

    /*
    
        lissage laplacien facteur de diffusion n iter
    
    */
    auto newMesh3 = geomAlgoLib::facteur_de_diffusion_n_iter(myMesh,1,0.5);

    geomAlgoLib::writeOFF(newMesh3,"output4.off");

    /*
    
        lissage laplacien de taubin n iter
    
    */

    auto newMesh4 = geomAlgoLib::lissage_de_taubin_n_iter(myMesh,50,0.33,-0.34);  //2.0 et -1.5 pour 100 iteration est intéressant

    geomAlgoLib::writeOFF(newMesh4,"output5.off");

    /*
    
        lissage laplacian ponderer via distances du vecteur v - u n iter
    
    */

    auto newMesh5 = geomAlgoLib::lissage_n_iter_ponderer(myMesh,10);

    geomAlgoLib::writeOFF(newMesh5,"output6.off");

    /*
    
        lissage laplacian ponderer via cotangentes n iter
    
    */

    auto newMesh6 = geomAlgoLib::lissage_n_iter_ponderer_cotangentes(myMesh,5);

    geomAlgoLib::writeOFF(newMesh6,"output7.off");

    std::cout << "The end..." << std::endl;

    return 0;
}