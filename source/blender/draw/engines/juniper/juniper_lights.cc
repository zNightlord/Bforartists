
#include "juniper_lights.hh"
#include "DNA_light_types.h"

namespace blender::juniper {

void Light::sync(const Object *ob) {
  const ::Light *la = (const ::Light *)ob->data;
  normalize_m4_m4(this->object_mat.ptr(), ob->object_to_world);

  this->light_groups = la->light_group_bits;
  if (la->type == LA_SUN) {
    this->type = LIGHT_TYPE_SUN;
  } else {
    this->type = LIGHT_TYPE_POINT;
  }

  this->color = float3(&la->r);
  this->power = la->energy;
  this->radius = la->radius;
}

void LightManager::begin_sync() {
  // Reset total light count
  meta_buf.num_lights = 0;
}

// TODO this should use DrawData handles, not base object pointer
// Since for shadows we will want the recalc flag to avoid updating every frame
void LightManager::sync_light(Object *ob) {
  Light &light = scene_lights.lookup_or_add_default(ob);
  light.used = true;
  meta_buf.num_lights++;

  // Sync light data with blender light
  light.sync(ob);
}

void LightManager::end_sync() {
  // Delete any unused lights
  Vector<Object*, 0> deleted_objects;

  // Resize to number of actually used lights
  light_buf.resize(meta_buf.num_lights);

  int light_idx = 0;
  for (auto item : scene_lights.items()) {
    Light &light = item.value;
    if (!light.used) {
      deleted_objects.append(item.key);
      continue;
    }

    // Add to buffer
    light_buf[light_idx++] = light;

    // Reset for next frame
    light.used = false;
  }

  // Update GPU buffers
  light_buf.push_update();
  meta_buf.push_update();

  for (auto &ob : deleted_objects) {
    scene_lights.remove(ob);
  }
}


}
