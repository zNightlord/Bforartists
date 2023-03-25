#pragma once
#include "draw_command.hh"
#include "GPU_batch.h"
#include "gpu_batch_private.hh"


namespace blender::bnpr
{
    struct GeomExtractionData
    {
        GPUBatch *batch;
        uint instance_len;
        uint vertex_len;
        uint vertex_first;
    };
    
    class GeomExtractCommand
    {
    public:
        /**
         * \brief see GPU_batch_draw_parameter_get and DrawCommandBuf::bind
         */
        static void gpu_batch_draw_parameter_get(
            GPUBatch *gpu_batch,
            int *batch_vert_len, int *batch_vert_first,
            int *batch_base_index,
            int *batch_inst_len
        );
    private:
        

        
    };
}
