#include <iostream>

#include "meshslicing.h"
#include "subdivloop.h"
#include "subdivlinear.h"
#include "elapsedtimer.h"
#include "oversampling.h"

REMESHER_NAMESPACE_BEGIN

void
Oversampling::start_oversampling (void)
{
  if (this->config.use_loop_subdiv)
  {
    this->run_loop_subdiv();
  }

  if (this->config.use_linear_subdiv)
  {
    this->run_linear_subdiv();
  }

  if (this->config.use_slicing)
  {
    this->run_grid_slicing();
  }
}

/* ---------------------------------------------------------------- */

void
Oversampling::run_loop_subdiv (void)
{
  SubdivLoop subdiv;
  subdiv.set_mesh(this->mesh);
  subdiv.set_vertex_info(this->vinfo);
  subdiv.set_feature_edges(this->features);
  subdiv.start_subdiv();
}

/* ---------------------------------------------------------------- */

void
Oversampling::run_linear_subdiv (void)
{
  /* Invalidate features if set. */
  if (this->features.get() != 0)
    this->features->clear();

  SubdivLinear subdiv;
  subdiv.set_mesh(this->mesh);
  subdiv.set_vertex_info(this->vinfo);
  subdiv.start_subdiv();
}

/* ---------------------------------------------------------------- */

void
Oversampling::run_grid_slicing (void)
{
  /* Invalidate features if set. */
  if (this->features.get() != 0)
    this->features->clear();

  /* Setup mesh slicer. */
  MeshSlicingConf msc;
  msc.vertex_snapping = 0.0001f;
  msc.edge_snapping = 0.2f;

  MeshSlicing ms;
  ms.set_config(msc);
  ms.set_mesh(this->mesh);

  /* Setup grid slicing. */
  Vec3f axis[3];
  axis[0] = Vec3f(1.0f, 0.0f, 0.0f);
  axis[1] = Vec3f(0.0f, 1.0f, 0.0f);
  axis[2] = Vec3f(0.0f, 0.0f, 1.0f);

  std::size_t slices[3];
  slices[0] = this->config.slices_x;
  slices[1] = this->config.slices_y;
  slices[2] = this->config.slices_z;

  for (std::size_t i = 0; i < 3; ++i)
  {
    Vec3f const& pn = axis[i];
    std::size_t num = slices[i];

    for (std::size_t j = 1; j <= num; ++j)
    {
      std::cout << "Slicing layer " << j << " of " << num << std::endl;
      float pos = (float)j / ((float)num + 1) - 0.5f;

      /* Determine if we use lazy or immediate triangulation. */
      if (this->config.lazy_triangulation)
        ms.lazy_slice(Plane3f(pn, pos));
      else
        ms.slice(Plane3f(pn, pos));
    }
  }

  /* Finally triangulate if lazy slicing is enabled. */
  if (this->config.lazy_triangulation)
    ms.triangulate();
}

REMESHER_NAMESPACE_END
