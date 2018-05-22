#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include "libremesh/modelwriter.h"
#include "libremesh/interface.h"

int main (int argc, char** argv) {

    // Parse command line arguments
    if (argc < 5) {
        std::cout << "Usage: " << argv[0] << " <input mesh> <output mesh> <mesh_verts> <lloyd_iters>" << std::endl;
        return 1;
    }
    std::string infile(argv[1]);
    std::string outfile(argv[2]);
    int mesh_verts = atoi(argv[3]);
    int lloyd_iters = atoi(argv[4]);

    // Loading, cleaning
    Remesher::Interface iface;
    {
        iface.load_model(infile);
        iface.clean_reference_mesh();
        iface.optimize_reference_mesh();
        Remesher::TriangleMeshPtr mesh = iface.get_reference_mesh();
    }

    // Resampling
    {
        Remesher::ResamplingConf resampling_conf;
        resampling_conf.sample_amount = mesh_verts;
        resampling_conf.perform_decimation = true;
        iface.set_resampling_conf(resampling_conf);
        iface.exec_resampling();
    }

    // Lloyd
    {
        Remesher::RelaxationConf lloyd_conf;
        lloyd_conf.iterations = lloyd_iters;
        iface.set_lloyd_conf(lloyd_conf);
        iface.exec_lloyd();
    }

    // Saving
    {
        Remesher::TriangleMeshPtr mesh(iface.get_evolving_mesh());
        Remesher::ModelWriter::save_model(outfile, mesh);
    }

    return 0;
}

