#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include "microdelaunay.h"
#include "cvtstats.h"

REMESHER_NAMESPACE_BEGIN

void
CvtStats::calc_stats (std::string const& gp_filename)
{
  MeshVertexList const& everts(this->emesh->get_vertices());

  //if (this->evinfo.get() == 0 || this->evinfo->size() != everts.size())
  //  this->evinfo = VertexInfoList::create(this->emesh);

  std::vector<std::pair<float, float> > meshinfo;

  /* Calc the mass for each Voronoi region. */
  for (std::size_t i = 0; i < everts.size(); ++i)
  {
    if (this->evinfo[i].vclass != VERTEX_CLASS_SIMPLE)
      continue;

    VertexInfo::VertexList vlist;
    this->evinfo->get_adjacent_vertices(i, vlist);

    /* Build micro patch. */
    MicroDelaunay md;
    md.set_center_vertex(everts[i]);
    for (std::size_t j = 0; j < vlist.size(); ++j)
      md.append_adjacent_vertex(everts[vlist[j]]);
    md.calculate_patch();

    /* Set density values if available. */
    if (this->rdens.get() != 0 && !this->rdens->empty())
    {
      bool valid_density = true;
      {
        VertexRef const& ref = this->evref[i];
        VertexDensity vd = this->rdens->get_density(ref.face, ref.bary);
        md.set_center_density(vd.density);
        valid_density = valid_density && vd.valid;
      }
      for (std::size_t i = 0; i < vlist.size(); ++i)
      {
        VertexRef const& ref = this->evref[vlist[i]];
        VertexDensity vd = this->rdens->get_density(ref.face, ref.bary);
        md.append_adjacent_density(vd.density);
        valid_density = valid_density && vd.valid;
      }
      if (!valid_density)
        md.clear_adjacent_density();
    }

    /* Get mass of the Voronoi cell. */
    float area, mass;
    md.get_voronoi_area_and_mass(area, mass);
    meshinfo.push_back(std::make_pair(area, mass));
  }

  if (gp_filename.empty())
    return;

  std::ofstream out_a((gp_filename + "_area").c_str());
  std::ofstream out_m((gp_filename + "_mass").c_str());
  for (std::size_t i = 0; i < meshinfo.size(); ++i)
  {
    out_a << i << " " << meshinfo[i].first << std::endl;
    out_m << i << " " << meshinfo[i].second << std::endl;
  }
  std::cout << "Wrote " << meshinfo.size() << " values to file" << std::endl;
  out_a.close();
  out_m.close();
}

REMESHER_NAMESPACE_END
