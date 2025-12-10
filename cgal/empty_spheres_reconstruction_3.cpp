#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Maximal_empty_spheres/contact_points_from_signed_distances.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <string>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef Kernel::Sphere_3 Sphere;
typedef std::pair<Point, Vector> Point_with_normal;

int main(int argc, char** argv) {

    bool filter_contact_spheres_with_bbox = true;
    int debug_level = 1;

    const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("data/3D/spheres.csv");;
    std::ifstream in(filename);

    CGAL::Timer timer;
    std::vector<std::pair<Sphere,int>> input_spheres;
    std::vector<Point_with_normal> pwns;

    double x, y, z, r;
    while(in >> x){
        in.ignore(10,','); in >> y;  in.ignore(10,','); in >> z; in.ignore(10,','); in >> r;
        int inside = r < 0 ? -1 : 1; // positive radius means outside, negative radius means inside
        input_spheres.emplace_back(std::make_pair(Sphere(Point(x,y,z),CGAL::square(r)),inside));
    }
    std::cout << "Read " << input_spheres.size() << " spheres" << std::endl;

    timer.start();
    CGAL::contact_points_from_signed_distances(input_spheres, std::back_inserter(pwns), filter_contact_spheres_with_bbox, debug_level);

    std::cout << "Computed " << pwns.size() << " contact points with normals in " << timer.time() << " sec." << std::endl;

    std::ofstream feil;
	feil.open("pwn.csv");
    for (int i=0; i<pwns.size(); i++){
        feil << pwns[i].first[0]  << "," << pwns[i].first[1]  << "," << pwns[i].first[2]  << ","
             << pwns[i].second[0] << "," << pwns[i].second[1] << "," << pwns[i].second[2] << std::endl;
    }
    feil.close();

    return 0;
}
