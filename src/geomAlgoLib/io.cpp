#include "io.hpp"

namespace geomAlgoLib
{

bool readOFF(const std::string& filePath, Mesh& mesh){

    std::ifstream input(filePath);

	if (!input || !(input >> mesh) || mesh.is_empty())
	{
		std::cerr << "Invalid .OFF file" << std::endl;
		return false;
	}

    return true;

}

void writeOFF(const Mesh& mesh, const std::string& filePath)
{
	std::ofstream in_myfile;
	in_myfile.open(filePath);

	CGAL::set_ascii_mode(in_myfile);

	in_myfile << "COFF" << std::endl // "COFF" makes the file support color information
			  << mesh.size_of_vertices() << ' ' 
			  << mesh.size_of_facets() << " 0" << std::endl; 
			  // nb of vertices, faces and edges (the latter is optional, thus 0)

	std::copy(mesh.points_begin(), mesh.points_end(),
			  std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));

	for (FacetCstIt i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
		HalfedgeFacetCstCirc j = i->facet_begin();

		CGAL_assertion(CGAL::circulator_size(j) >= 3);

		in_myfile << CGAL::circulator_size(j) << ' ';
		do
		{
			in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());

		} while (++j != i->facet_begin());

		in_myfile << std::endl;
	}

	in_myfile.close();

	std::cout << "Successfully exported at path: " << filePath << " !" << std::endl;
}

void StoreinFiles(const Mesh& mesh,const FacetDoubleMap & tab, const std::string& filePath ){
	std::ofstream in_myfile;
	in_myfile.open(filePath);

	CGAL::set_ascii_mode(in_myfile);

	in_myfile << "COFF" << std::endl // "COFF" makes the file support color information
			  << mesh.size_of_vertices() << ' ' 
			  << mesh.size_of_facets() << ' ' << mesh.size_of_halfedges() /2 << std::endl; 
			  // nb of vertices, faces and edges (the latter is optional, thus 0)

	std::copy(mesh.points_begin(), mesh.points_end(),
			  std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));


	for (auto i = tab.begin(); i != tab.end(); ++i)
	{
		HalfedgeFacetCstCirc j = i->first->facet_begin();

		CGAL_assertion(CGAL::circulator_size(j) >= 3);

		in_myfile << CGAL::circulator_size(j) << ' ';
		do
		{
			in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());
			

		} while (++j != i->first->facet_begin());
		do
		{
			int c = geomAlgoLib::calcul_couleur_mesh(i->second);
			in_myfile << ' ' << c  ;
			

		} while (++j != i->first->facet_begin());
		in_myfile << std::endl;
	}

	in_myfile.close();

	std::cout << "Successfully exported at path: " << filePath << " !" << std::endl;
}


double calcul_couleur_mesh(const double& second){
	int c = 255 / (second + 1);   
	return c;
   }

void StoreinFiles_labels(const Mesh& mesh,const FacetStringMap & tab1,const FacetDoubleMap & tab2 ,const std::string& filePath ){
	std::ofstream in_myfile;
	in_myfile.open(filePath);

	CGAL::set_ascii_mode(in_myfile);

	in_myfile << "COFF" << std::endl // "COFF" makes the file support color information
			  << mesh.size_of_vertices() << ' ' 
			  << mesh.size_of_facets() << ' ' << mesh.size_of_halfedges() /2 << std::endl; 
			  // nb of vertices, faces and edges (the latter is optional, thus 0)

	std::copy(mesh.points_begin(), mesh.points_end(),
			  std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));


	for (auto i = tab1.begin(); i != tab1.end(); ++i)
	{
		HalfedgeFacetCstCirc j = i->first->facet_begin();

		CGAL_assertion(CGAL::circulator_size(j) >= 3);

		in_myfile << CGAL::circulator_size(j) << ' ';
		do
		{
			in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());
			

		} 
		while (++j != i->first->facet_begin());
		if (i->second == "Grande Face"){
			in_myfile << ' ' << 255 << ' ' << 255 << ' ' << 255  ;
			
		}
		else if(i->second == "orientées vers le haut"){
			in_myfile << ' ' << 0 << ' ' << 0 << ' ' << 255;
		}	
		else if(i->second == "obstacle"){
			in_myfile << ' ' << 255 << ' ' << 0 << ' ' << 0;
		}
		else if(i->second == "a portee"){
			in_myfile << ' ' << 0 << ' ' << 255 << ' ' << 0;
		}
		else{
		do
		{
		int c = geomAlgoLib::calcul_couleur_mesh(tab2.at(i->first));
			
			
		in_myfile << ' ' << c  ;
			

		} while (++j != i->first->facet_begin());
			
		}
			in_myfile << std::endl;
			
		}
		
		
	

	in_myfile.close();

	std::cout << "Successfully exported at path: " << filePath << " !" << std::endl;
}

}