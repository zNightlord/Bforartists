/* SPDX-FileCopyrightText: 2023 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup RNA
 */

#include <cstdlib>

#include "RNA_define.hh"
#include "RNA_enum_types.hh"

#include "rna_internal.hh"

#include "DNA_lightprobe_types.h"

#include "WM_types.hh"

#ifdef RNA_RUNTIME

#  include "MEM_guardedalloc.h"

#  include "BKE_main.hh"
#  include "DEG_depsgraph.hh"

#  include "DNA_collection_types.h"
#  include "DNA_object_types.h"

#  include "WM_api.hh"

static StructRNA *rna_LightProbe_refine(PointerRNA *ptr)
{
  LightProbe *probe = (LightProbe *)ptr->data;
  switch (probe->type) {
    case LIGHTPROBE_TYPE_PLANE:
      return &RNA_LightProbePlane;
    case LIGHTPROBE_TYPE_SPHERE:
      return &RNA_LightProbeSphere;
    case LIGHTPROBE_TYPE_VOLUME:
      return &RNA_LightProbeVolume;
    default:
      return &RNA_LightProbe;
  }
}

static void rna_LightProbe_recalc(Main * /*bmain*/, Scene * /*scene*/, PointerRNA *ptr)
{
  DEG_id_tag_update(ptr->owner_id, ID_RECALC_GEOMETRY);
}

#else

static EnumPropertyItem parallax_type_items[] = {
    {LIGHTPROBE_SHAPE_ELIPSOID, "ELIPSOID", ICON_NONE, "Sphere", ""},
    {LIGHTPROBE_SHAPE_BOX, "BOX", ICON_NONE, "Box", ""},
    {0, nullptr, 0, nullptr, nullptr},
};

static EnumPropertyItem lightprobe_type_items[] = {
    {LIGHTPROBE_TYPE_SPHERE,
     "SPHERE",
     ICON_LIGHTPROBE_SPHERE,
     "Sphere",
     "Light probe that captures precise lighting from all directions at a single point in space"},
    {LIGHTPROBE_TYPE_PLANE,
     "PLANE",
     ICON_LIGHTPROBE_PLANE,
     "Plane",
     "Light probe that captures incoming light from a single direction on a plane"},
    {LIGHTPROBE_TYPE_VOLUME,
     "VOLUME",
     ICON_LIGHTPROBE_VOLUME,
     "Volume",
     "Light probe that captures low frequency lighting inside a volume"},
    {0, nullptr, 0, nullptr, nullptr},
};

static void rna_def_lightprobe(BlenderRNA *brna)
{
  StructRNA *srna;
  PropertyRNA *prop;

  srna = RNA_def_struct(brna, "LightProbe", "ID");
  RNA_def_struct_refine_func(srna, "rna_LightProbe_refine");
  RNA_def_struct_ui_text(
      srna, "LightProbe", "Light Probe data for lighting capture objects"); /*BFA - not data-block*/
  RNA_def_struct_ui_icon(srna, ICON_OUTLINER_DATA_LIGHTPROBE);

  prop = RNA_def_property(srna, "type", PROP_ENUM, PROP_NONE);
  RNA_def_property_enum_items(prop, lightprobe_type_items);
  RNA_def_property_ui_text(prop, "Type", "Type of light probe");
  RNA_def_property_clear_flag(prop, PROP_EDITABLE);

  prop = RNA_def_property(srna, "clip_start", PROP_FLOAT, PROP_DISTANCE);
  RNA_def_property_float_sdna(prop, nullptr, "clipsta");
  RNA_def_property_range(prop, 1e-6f, FLT_MAX);
  RNA_def_property_ui_range(prop, 0.001f, FLT_MAX, 10, 3);
  RNA_def_property_ui_text(
      prop, "Clip Start", "Probe clip start, below which objects will not appear in reflections");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  prop = RNA_def_property(srna, "show_clip", PROP_BOOLEAN, PROP_NONE);
  RNA_def_property_boolean_sdna(prop, nullptr, "flag", LIGHTPROBE_FLAG_SHOW_CLIP_DIST);
  RNA_def_property_ui_text(prop, "Clipping", "Show the clipping distances in the 3D view");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, nullptr);

  prop = RNA_def_property(srna, "show_influence", PROP_BOOLEAN, PROP_NONE);
  RNA_def_property_boolean_sdna(prop, nullptr, "flag", LIGHTPROBE_FLAG_SHOW_INFLUENCE);
  RNA_def_property_ui_text(prop, "Influence", "Show the influence volume in the 3D view");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, nullptr);

  prop = RNA_def_property(srna, "influence_distance", PROP_FLOAT, PROP_DISTANCE);
  RNA_def_property_float_sdna(prop, nullptr, "distinf");
  RNA_def_property_range(prop, 0.0f, FLT_MAX);
  RNA_def_property_ui_text(prop, "Influence Distance", "Influence distance of the probe");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, nullptr);

#  if 1 /* Deprecated: Remove in Blender 4.5 */
  prop = RNA_def_property(srna, "visibility_buffer_bias", PROP_FLOAT, PROP_NONE);
  RNA_def_property_float_sdna(prop, nullptr, "vis_bias");
  RNA_def_property_range(prop, 0.001f, 9999.0f);
  RNA_def_property_ui_range(prop, 0.001f, 5.0f, 1.0, 3);
  RNA_def_property_ui_text(
      prop, "Visibility Bias", "Bias for reducing self shadowing (Deprecated)");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, nullptr);

  prop = RNA_def_property(srna, "visibility_bleed_bias", PROP_FLOAT, PROP_FACTOR);
  RNA_def_property_float_sdna(prop, nullptr, "vis_bleedbias");
  RNA_def_property_range(prop, 0.0f, 1.0f);
  RNA_def_property_ui_text(prop,
                           "Visibility Bleed Bias",
                           "Bias for reducing light-bleed on variance shadow maps (Deprecated)");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, nullptr);

  prop = RNA_def_property(srna, "visibility_blur", PROP_FLOAT, PROP_FACTOR);
  RNA_def_property_float_sdna(prop, nullptr, "vis_blur");
  RNA_def_property_range(prop, 0.0f, 1.0f);
  RNA_def_property_ui_text(
      prop, "Visibility Blur", "Filter size of the visibility blur (Deprecated)");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  prop = RNA_def_property(srna, "visibility_collection", PROP_POINTER, PROP_NONE);
  RNA_def_property_struct_type(prop, "Collection");
  RNA_def_property_pointer_sdna(prop, nullptr, "visibility_grp");
  RNA_def_property_flag(prop, PROP_EDITABLE);
  RNA_def_property_override_flag(prop, PROPOVERRIDE_OVERRIDABLE_LIBRARY);
  RNA_def_property_ui_text(
      prop, "Visibility Collection", "Restrict objects visible for this probe (Deprecated)");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  prop = RNA_def_property(srna, "invert_visibility_collection", PROP_BOOLEAN, PROP_NONE);
  RNA_def_property_boolean_sdna(prop, nullptr, "flag", LIGHTPROBE_FLAG_INVERT_GROUP);
  RNA_def_property_flag(prop, PROP_EDITABLE);
  RNA_def_property_ui_text(prop, "Invert Collection", "Invert visibility collection (Deprecated)");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");
#  endif

  /* Data preview */
  prop = RNA_def_property(srna, "show_data", PROP_BOOLEAN, PROP_NONE);
  RNA_def_property_boolean_sdna(prop, nullptr, "flag", LIGHTPROBE_FLAG_SHOW_DATA);
  RNA_def_property_ui_icon(prop, ICON_HIDE_ON, 1);
  RNA_def_property_ui_text(
      prop, "Display Data (Deprecated)", "Deprecated, use use_data_display instead");
  RNA_def_property_override_flag(prop, PROPOVERRIDE_OVERRIDABLE_LIBRARY);
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, nullptr);

  prop = RNA_def_property(srna, "use_data_display", PROP_BOOLEAN, PROP_NONE);
  RNA_def_property_boolean_sdna(prop, nullptr, "flag", LIGHTPROBE_FLAG_SHOW_DATA);
  RNA_def_property_ui_text(
      prop, "Display Data", "Display sampled data in the viewport to debug captured light");
  RNA_def_property_override_flag(prop, PROPOVERRIDE_OVERRIDABLE_LIBRARY);
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, nullptr);

  prop = RNA_def_property(srna, "data_display_size", PROP_FLOAT, PROP_FACTOR);
  RNA_def_property_float_sdna(prop, nullptr, "data_display_size");
  RNA_def_property_range(prop, 0.0f, FLT_MAX);
  RNA_def_property_ui_range(prop, 0.01f, 1.0f, 1, 3);
  RNA_def_property_ui_text(prop, "Display Data Size", "Viewport display size of the sampled data");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, nullptr);

  /* common */
  rna_def_animdata_common(srna);
}

static void rna_def_lightprobe_plane(BlenderRNA *brna)
{
  StructRNA *srna;

  srna = RNA_def_struct(brna, "LightProbePlane", "LightProbe");
  RNA_def_struct_sdna(srna, "LightProbe");
  RNA_def_struct_ui_text(
      srna,
      "Planar Probe",
      "Light probe that captures incoming light from a single direction on a plane");
  RNA_def_struct_ui_icon(srna, ICON_LIGHTPROBE_PLANE);
}

static void rna_def_lightprobe_sphere(BlenderRNA *brna)
{
  StructRNA *srna;
  PropertyRNA *prop;

  srna = RNA_def_struct(brna, "LightProbeSphere", "LightProbe");
  RNA_def_struct_sdna(srna, "LightProbe");
  RNA_def_struct_ui_text(
      srna,
      "Spherical Probe",
      "Light probe that captures precise lighting from all directions at a single point in space");
  RNA_def_struct_ui_icon(srna, ICON_LIGHTPROBE_SPHERE);

  prop = RNA_def_property(srna, "influence_type", PROP_ENUM, PROP_NONE);
  RNA_def_property_enum_sdna(prop, nullptr, "attenuation_type");
  RNA_def_property_enum_items(prop, parallax_type_items);
  RNA_def_property_ui_text(prop, "Type", "Type of influence volume");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, nullptr);

  prop = RNA_def_property(srna, "falloff", PROP_FLOAT, PROP_FACTOR);
  RNA_def_property_range(prop, 0.0f, 1.0f);
  RNA_def_property_ui_text(prop, "Falloff", "Control how fast the probe influence decreases");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, nullptr);

  prop = RNA_def_property(srna, "clip_end", PROP_FLOAT, PROP_DISTANCE);
  RNA_def_property_float_sdna(prop, nullptr, "clipend");
  RNA_def_property_range(prop, 1e-6f, FLT_MAX);
  RNA_def_property_ui_range(prop, 0.001f, FLT_MAX, 10, 3);
  RNA_def_property_ui_text(
      prop, "Clip End", "Probe clip end, beyond which objects will not appear in reflections");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  /* Custom parallax */
  prop = RNA_def_property(srna, "use_custom_parallax", PROP_BOOLEAN, PROP_NONE);
  RNA_def_property_boolean_sdna(prop, nullptr, "flag", LIGHTPROBE_FLAG_CUSTOM_PARALLAX);
  RNA_def_property_ui_text(
      prop, "Use Custom Parallax", "Enable custom settings for the parallax correction volume");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, nullptr);

  prop = RNA_def_property(srna, "show_parallax", PROP_BOOLEAN, PROP_NONE);
  RNA_def_property_boolean_sdna(prop, nullptr, "flag", LIGHTPROBE_FLAG_SHOW_PARALLAX);
  RNA_def_property_ui_text(prop, "Parallax", "Show the parallax correction volume in the 3D view\nThe Parallax checkbox turns active when Custom Parallax is activated"); /* BFA */
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, nullptr);

  prop = RNA_def_property(srna, "parallax_type", PROP_ENUM, PROP_NONE);
  RNA_def_property_enum_items(prop, parallax_type_items);
  RNA_def_property_ui_text(prop, "Type", "Type of parallax volume");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, nullptr);

  prop = RNA_def_property(srna, "parallax_distance", PROP_FLOAT, PROP_DISTANCE);
  RNA_def_property_float_sdna(prop, nullptr, "distpar");
  RNA_def_property_range(prop, 0.0f, FLT_MAX);
  RNA_def_property_ui_text(prop, "Parallax Radius", "Lowest corner of the parallax bounding box");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, nullptr);
}

static void rna_def_lightprobe_volume(BlenderRNA *brna)
{
  StructRNA *srna;
  PropertyRNA *prop;

  srna = RNA_def_struct(brna, "LightProbeVolume", "LightProbe");
  RNA_def_struct_sdna(srna, "LightProbe");
  RNA_def_struct_ui_text(
      srna, "Volume Probe", "Light probe that captures low frequency lighting inside a volume");
  RNA_def_struct_ui_icon(srna, ICON_LIGHTPROBE_VOLUME);

  prop = RNA_def_property(srna, "intensity", PROP_FLOAT, PROP_NONE);
  RNA_def_property_float_sdna(prop, nullptr, "intensity");
  RNA_def_property_range(prop, 0.0f, FLT_MAX);
  RNA_def_property_ui_range(prop, 0.0f, 3.0f, 1.0, 3);
  RNA_def_property_ui_text(
      prop, "Intensity", "Modify the intensity of the lighting captured by this probe");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  /* Irradiance grid */
  prop = RNA_def_property(srna, "resolution_x", PROP_INT, PROP_NONE);
  RNA_def_property_int_sdna(prop, nullptr, "grid_resolution_x");
  RNA_def_property_range(prop, 1, 256);
  RNA_def_property_ui_text(
      prop, "Resolution X", "Number of samples along the x axis of the volume");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  prop = RNA_def_property(srna, "resolution_y", PROP_INT, PROP_NONE);
  RNA_def_property_int_sdna(prop, nullptr, "grid_resolution_y");
  RNA_def_property_range(prop, 1, 256);
  RNA_def_property_ui_text(
      prop, "Resolution Y", "Number of samples along the y axis of the volume");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  prop = RNA_def_property(srna, "resolution_z", PROP_INT, PROP_NONE);
  RNA_def_property_int_sdna(prop, nullptr, "grid_resolution_z");
  RNA_def_property_range(prop, 1, 256);
  RNA_def_property_ui_text(
      prop, "Resolution Z", "Number of samples along the z axis of the volume");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  /* NOTE: We reuse the same DNA for this property for historical reason, but we want a different
   * name and tooltip for it. */
  prop = RNA_def_property(srna, "capture_distance", PROP_FLOAT, PROP_DISTANCE);
  RNA_def_property_float_sdna(prop, nullptr, "clipend");
  RNA_def_property_range(prop, 1e-6f, FLT_MAX);
  RNA_def_property_ui_range(prop, 0.001f, FLT_MAX, 10, 1);
  RNA_def_property_ui_text(prop,
                           "Capture Distance",
                           "Distance around the probe volume that will be considered "
                           "during the bake");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  prop = RNA_def_property(srna, "normal_bias", PROP_FLOAT, PROP_FACTOR);
  RNA_def_property_float_sdna(prop, nullptr, "grid_normal_bias");
  RNA_def_property_ui_text(prop,
                           "Normal Bias",
                           "Offset sampling of the irradiance grid in "
                           "the surface normal direction to reduce light bleeding");
  RNA_def_property_range(prop, 0.0f, FLT_MAX);
  RNA_def_property_ui_range(prop, 0.0f, 1.0f, 1, 2);
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  prop = RNA_def_property(srna, "view_bias", PROP_FLOAT, PROP_FACTOR);
  RNA_def_property_float_sdna(prop, nullptr, "grid_view_bias");
  RNA_def_property_ui_text(prop,
                           "View Bias",
                           "Offset sampling of the irradiance grid in "
                           "the viewing direction to reduce light bleeding");
  RNA_def_property_range(prop, 0.0f, FLT_MAX);
  RNA_def_property_ui_range(prop, 0.0f, 1.0f, 1, 2);
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  prop = RNA_def_property(srna, "facing_bias", PROP_FLOAT, PROP_FACTOR);
  RNA_def_property_float_sdna(prop, nullptr, "grid_facing_bias");
  RNA_def_property_ui_text(
      prop, "Facing Bias", "Smoother irradiance interpolation but introduce light bleeding");
  RNA_def_property_range(prop, 0.0f, FLT_MAX);
  RNA_def_property_ui_range(prop, 0.0f, 1.0f, 1, 2);
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  prop = RNA_def_property(srna, "bake_samples", PROP_INT, PROP_NONE);
  RNA_def_property_int_sdna(prop, nullptr, "grid_bake_samples");
  RNA_def_property_ui_text(
      prop, "Bake Samples", "Number of ray directions to evaluate when baking");
  RNA_def_property_range(prop, 1, INT_MAX);
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  prop = RNA_def_property(srna, "surface_bias", PROP_FLOAT, PROP_FACTOR);
  RNA_def_property_float_sdna(prop, nullptr, "grid_surface_bias");
  RNA_def_property_ui_text(
      prop, "Surface Offset", "Moves capture points away from surfaces to prevent artifacts");
  RNA_def_property_range(prop, 0.0f, 1.0f);
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  prop = RNA_def_property(srna, "escape_bias", PROP_FLOAT, PROP_FACTOR);
  RNA_def_property_float_sdna(prop, nullptr, "grid_escape_bias");
  RNA_def_property_ui_text(prop,
                           "Search Distance",
                           "Distance to search for valid capture positions to prevent "
                           "lighting artifacts");
  RNA_def_property_range(prop, 0.0f, 1.0f);
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  prop = RNA_def_property(srna, "surfel_density", PROP_INT, PROP_NONE);
  RNA_def_property_int_sdna(prop, nullptr, "grid_surfel_density");
  RNA_def_property_range(prop, 1, INT_MAX);
  RNA_def_property_ui_text(
      prop,
      "Surfel Resolution",
      "Number of surfels to spawn in one local unit distance (higher values improve quality)");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  prop = RNA_def_property(srna, "validity_threshold", PROP_FLOAT, PROP_FACTOR);
  RNA_def_property_float_sdna(prop, nullptr, "grid_validity_threshold");
  RNA_def_property_ui_text(prop,
                           "Validity Threshold",
                           "Ratio of front-facing surface hits under which a grid sample will "
                           "not be considered for lighting");
  RNA_def_property_range(prop, 0.0f, 1.0f);
  RNA_def_property_ui_range(prop, 0.0f, 1.0f, 1, 2);
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  prop = RNA_def_property(srna, "dilation_threshold", PROP_FLOAT, PROP_FACTOR);
  RNA_def_property_float_sdna(prop, nullptr, "grid_dilation_threshold");
  RNA_def_property_ui_text(prop,
                           "Dilation Threshold",
                           "Ratio of front-facing surface hits under which a grid sample will "
                           "reuse neighbors grid sample lighting");
  RNA_def_property_range(prop, 0.0f, 1.0f);
  RNA_def_property_ui_range(prop, 0.0f, 1.0f, 1, 2);
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  prop = RNA_def_property(srna, "dilation_radius", PROP_FLOAT, PROP_FACTOR);
  RNA_def_property_float_sdna(prop, nullptr, "grid_dilation_radius");
  RNA_def_property_ui_text(
      prop,
      "Dilation Radius",
      "Radius in grid sample to search valid grid samples to copy into invalid grid samples");
  RNA_def_property_range(prop, 1.0f, 5.0f);
  RNA_def_property_ui_range(prop, 1.0f, 5.0f, 1, 2);
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  prop = RNA_def_property(srna, "capture_world", PROP_BOOLEAN, PROP_NONE);
  RNA_def_property_boolean_sdna(prop, nullptr, "grid_flag", LIGHTPROBE_GRID_CAPTURE_WORLD);
  RNA_def_property_ui_text(
      prop,
      "Capture World",
      "Bake incoming light from the world instead of just the visibility "
      "for more accurate lighting, but lose correct blending to surrounding irradiance volumes");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, nullptr);

  prop = RNA_def_property(srna, "capture_indirect", PROP_BOOLEAN, PROP_NONE);
  RNA_def_property_boolean_sdna(prop, nullptr, "grid_flag", LIGHTPROBE_GRID_CAPTURE_INDIRECT);
  RNA_def_property_ui_text(prop,
                           "Capture Indirect",
                           "Bake light bounces from light sources for more accurate lighting");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, nullptr);

  prop = RNA_def_property(srna, "capture_emission", PROP_BOOLEAN, PROP_NONE);
  RNA_def_property_boolean_sdna(prop, nullptr, "grid_flag", LIGHTPROBE_GRID_CAPTURE_EMISSION);
  RNA_def_property_ui_text(
      prop, "Capture Emission", "Bake emissive surfaces for more accurate lighting");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, nullptr);

  prop = RNA_def_property(srna, "clamp_direct", PROP_FLOAT, PROP_NONE);
  RNA_def_property_float_sdna(prop, nullptr, "grid_clamp_direct");
  RNA_def_property_range(prop, 0.0f, FLT_MAX);
  RNA_def_property_ui_range(prop, 0.0f, FLT_MAX, 1, 2);
  RNA_def_property_ui_text(
      prop, "Clamp Direct", "Clamp the direct lighting intensity to reduce noise (0 to disable)");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");

  prop = RNA_def_property(srna, "clamp_indirect", PROP_FLOAT, PROP_NONE);
  RNA_def_property_float_sdna(prop, nullptr, "grid_clamp_indirect");
  RNA_def_property_range(prop, 0.0f, FLT_MAX);
  RNA_def_property_ui_range(prop, 0.0f, FLT_MAX, 1, 2);
  RNA_def_property_ui_text(prop,
                           "Clamp Indirect",
                           "Clamp the indirect lighting intensity to reduce noise (0 to disable)");
  RNA_def_property_update(prop, NC_MATERIAL | ND_SHADING, "rna_LightProbe_recalc");
}

void RNA_def_lightprobe(BlenderRNA *brna)
{
  rna_def_lightprobe(brna);
  rna_def_lightprobe_plane(brna);
  rna_def_lightprobe_sphere(brna);
  rna_def_lightprobe_volume(brna);
}

#endif
