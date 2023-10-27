#ifndef BNPR_COMPACTION__INCLUDED
#define BNPR_COMPACTION__INCLUDED

#if !defined(COMPACTION_LIB_EXCLUDE)
/* Input: */
/* #define COMPACTION_LIB_EXCLUDE --- define to remove this lib. used for separate kernels in the same file */
/* #define GLOBAL_COUNTER XXX --- Global Compaction Counter, should be cleared to 0*/
/* #define CP_TAG XXX --- Tag */

/* Macro expansion, for details, see
/* ---------------------------------------
/* https://stackoverflow.com/questions/1489932/how-to-concatenate-twice-with-the-c-preprocessor-and-expand-a-macro-as-in-arg */
#define CAT(x, y) CAT_(x, y)
#define CAT_(x, y) x ## y
 

#define LDS_DIGIT_PER_LANE CAT(LDS_digit_per_lane_, CP_TAG) 
shared uint LDS_DIGIT_PER_LANE[32u]; 

#define LDS_OFFSET_PER_LANE_SLOT CAT(LDS_offset_per_lane_slot_, CP_TAG)
shared uint LDS_OFFSET_PER_LANE_SLOT[32u];  

#define LDS_HIST_BLK CAT(LDS_hist_blk_, CP_TAG)
shared uint LDS_HIST_BLK;  

#define LDS_SCAN_BLOCK_OFFSET CAT(LDS_scan_block_offset_, CP_TAG)
shared uint LDS_SCAN_BLOCK_OFFSET; 


uint CAT(compact_, CP_TAG)(bool val, uint groupIdx)
{
	const uint wave_id = groupIdx >> 5u; /* must be < 32 which is ensured since tg size <= 1024 */
    const uint lane_id = groupIdx % 32u;
    const uint num_waves = gl_WorkGroupSize.x >> 5u; 
    
	/* Clear LDS counters */
    if (wave_id == 0u) 
    { 
        if (lane_id == 0u) LDS_HIST_BLK = 0u; 
        LDS_DIGIT_PER_LANE[lane_id] = 0u; 
    }
    barrier(); 
	/*  w0    w1    w2    w3         LDS_DIGIT_PER_LANE
	 * 00:1  04:1  08:0  12:1  l0    0000
	 * 01:0  05:0  09:0  13:0  l1    0000
	 * 02:1  06:1  10:0  14:0  l2    0000
	 * 03:1  07:0  11:1  15:0  l3    0000
	 * ---------------------------
	 * LDS_HIST_BLK
	 *  0
	*/

    /* Mark 1/0 at bit #wave_id */ 
    uint compact_bitval = uint(val); 
    uint lds_compact_input = compact_bitval << wave_id; 
    atomicOr(LDS_DIGIT_PER_LANE[lane_id], lds_compact_input); 
    barrier(); 
	/*  w0    w1    w2    w3         LDS_DIGIT_PER_LANE
	 * 00:1  04:1  08:0  12:1  l0    1011
	 * 01:0  05:0  09:0  13:0  l1    0000
	 * 02:1  06:1  10:0  14:0  l2    0011
	 * 03:1  07:0  11:1  15:0  l3    0101
	*/

    /* Prefix sum on lane sums */
    uint lane_digit = LDS_DIGIT_PER_LANE[lane_id]; 
    uint wave_mask = (~(0xffffffffu << wave_id));
	/*  w0    w1    w2    w3     
	 * 0000  0001  0011  0111        wave_mask
	 *
	 *  w0    w1    w2    w3         lane_digit
	 * 00:1  04:1  08:0  12:1  l0    1011
	 * 01:0  05:0  09:0  13:0  l1    0000
	 * 02:1  06:1  10:0  14:0  l2    0011
	 * 03:1  07:0  11:1  15:0  l3    0101
	*/
    uint lane_digit_masked = lane_digit & wave_mask; 
    uint num_1_bits_low = bitCount(lane_digit_masked);
    uint lane_offset = num_1_bits_low; 

    if (wave_id == 0u)
    {
        uint lane_digit_sum = bitCount(lane_digit); /* digit sum */ 
        LDS_OFFSET_PER_LANE_SLOT[lane_id] = atomicAdd(LDS_HIST_BLK, lane_digit_sum);  
    }
    barrier(); 

    /* Add block sum to global counter. */
    if (groupIdx == gl_WorkGroupSize.x - 1u)
    {
        LDS_SCAN_BLOCK_OFFSET = atomicAdd(
            GLOBAL_COUNTER, 
            LDS_HIST_BLK
        ); 
    }
    barrier(); 

    /* Compute final offset */
    uint local_offset = LDS_OFFSET_PER_LANE_SLOT[lane_id] + lane_offset; 
    uint blk_offset   = LDS_SCAN_BLOCK_OFFSET; 
    
    uint scanres = local_offset + blk_offset; 
	return scanres; 
}

#endif

#endif




