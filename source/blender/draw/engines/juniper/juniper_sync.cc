
#include "juniper_sync.hh"
#include "juniper_engine.h"
#include "BKE_object.h"
#include "draw_manager.hh"
#include "juniper_instance.hh"

namespace blender::juniper {

void SceneSync::sync_mesh(Object *ob) {
    draw::ObjectRef ref = DRW_object_ref_get(ob);
    draw::ResourceHandle handle = jnpr_.manager->resource_handle(ref);

    Vector<NPRMaterial *> &materials = jnpr_.materials.object_materials_get(ob);

    // Attrib usage can only be registered from one gpumat per surface segment
    // TODO find a solution for attribs used only in prepass
    Vector<GPUMaterial *> gpumats = Vector<GPUMaterial *>();
    for (NPRMaterial * mat : materials) {
        gpumats.append(mat->pass_main.gpumat);
        // gpumats.append(mat->pass_pre.gpumat);
        // gpumats.append(mat->pass_depth.gpumat);
    }

    GPUBatch **mat_geom = DRW_cache_object_surface_material_get(
            ob, gpumats.data(), gpumats.size());

    if (mat_geom == nullptr) {
        return;
    }

    for (auto i : materials.index_range()) {
        GPUBatch *geom = mat_geom[i];
        if (geom == nullptr) {
            continue;
        }

        NPRMaterial *material = materials[i];
        if (material->pass_pre.sub_pass != nullptr) {
          material->pass_pre.sub_pass->draw(geom, handle);
        }
        if (material->pass_main.sub_pass != nullptr) {
            material->pass_main.sub_pass->draw(geom, handle);
        }
    }

    // TODO add object attributes info
    jnpr_.manager->extract_object_attributes(handle, ref, gpumats);
}

void SceneSync::sync_light(Object *ob) {
  jnpr_.lights.sync_light(ob);
}

/** Add object to render*/
void SceneSync::sync_object(Object *ob) {
    // TODO support other object types
    const bool is_renderable_type = ELEM(ob->type, OB_MESH, OB_LAMP);

    const int ob_vis_flag = DRW_object_visibility_in_active_context(ob);

    const bool ob_is_visible = DRW_object_is_renderable(ob)
            && (ob_vis_flag & OB_VISIBLE_SELF) != 0;

    /* Unsupported or non-visible object. */
    if (!is_renderable_type || !ob_is_visible) {
        return;
    }

    switch (ob->type) {
        case OB_MESH:
            sync_mesh(ob);
            break;
        case OB_LAMP:
            sync_light(ob);
            break;
        default:
            return;
    }
}
}
