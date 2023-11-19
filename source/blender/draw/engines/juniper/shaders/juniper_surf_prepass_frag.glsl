
#pragma BLENDER_REQUIRE(juniper_nodetree_lib.glsl)

#define Closure float
float out_surface, out_volume, out_displacement, out_thickness;

// Forward declaration
vec4 nodetree_jnpr_color();

void main() {

    init_globals();
    // TODO normal from nodetree (?)
    float alpha = nodetree_jnpr_color().a;
    out_normal = safe_normalize(interp.N.xyz);

    if (alpha <= 0.0) {
        discard;
    }
}
