
#pragma BLENDER_REQUIRE(juniper_nodetree_lib.glsl)

// Forward declaration
vec4 nodetree_jnpr_color();

// Temp stubs since we're hack-using eevee output node.
// Should replace the node with something else lol
#define Closure float
float out_surface, out_volume, out_displacement, out_thickness;

void main() {
    init_globals();
    vec4 result = nodetree_jnpr_color();
    out_color.rgb = result.rgb;

    out_color.a = 1.0;
    if (result.a <= 0.0) {
        discard;
    }
}
