#include "bnpr_strokegen_geom_extract.hh"




// IB indexing:
// batch_base_index + batch_vert_first + vert index
// batch_vert_first + gl_VertexID
/**
 * \brief Get IB indexing info from a given gpu-batch.
 * \param gpu_batch gpu batch with mesh ib(s) & vb
 * \param batch_vert_len #verts in current batch
 * \param batch_vert_first ib offset within the batch's subrange
 * \param batch_base_index ib offset of verts of the batch,
 * need to add this with <cref> batch_vert_first unless we are read the index with gl_VertexID.
 * \param batch_inst_len #instances in current batch
 */
void blender::bnpr::GeomExtractCommand::gpu_batch_draw_parameter_get(
    GPUBatch* gpu_batch,
    int* batch_vert_len,
    int* batch_vert_first,
    int* batch_base_index,
    int* batch_inst_len
)
{
    const gpu::Batch* batch = static_cast<gpu::Batch*>(gpu_batch);

    if (batch->elem)
    {
        *batch_vert_len = batch->elem_()->index_len_get();
        *batch_vert_first = batch->elem_()->index_start_get();
        *batch_base_index = batch->elem_()->index_base_get();
    }
    else
    {
        *batch_vert_len = batch->verts_(0)->vertex_len;
        *batch_vert_first = 0;
        *batch_base_index = -1;
    }

    int i_count = (batch->inst[0]) ? batch->inst_(0)->vertex_len : 1;
    /* Meh. This is to be able to use different numbers of verts in instance VBO's. */
    if (batch->inst[1] != nullptr) { i_count = min_ii(i_count, batch->inst_(1)->vertex_len); }
    *batch_inst_len = i_count;
}
