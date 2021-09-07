// ************************************************************************
//
// miniAMR: stencil computations with boundary exchange and AMR.
//
// Copyright (2014) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
// retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
// Questions? Contact Courtenay T. Vaughan (ctvaugh@sandia.gov)
//                    Richard F. Barrett (rfbarre@sandia.gov)
//
// ************************************************************************

#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "block.h"
#include "comm.h"
#include "proto.h"
#include "timer.h"

// Profiling output.
void profile(void)
{
   int i;
   double total_gflops, gflops_rank, total_fp_ops, total_fp_adds,
          total_fp_divs, delta;
   object *op;
   char *version = "1.7.0";
   FILE *fp;

   calculate_results();
   total_fp_ops = average[139] + average[140] + average[141];
   total_gflops = total_fp_ops/(average[38]*1024.0*1024.0*1024.0);
   gflops_rank = total_gflops/((double) num_pes);

   if (!my_pe) {
      if (report_perf & 1) {
         fp = fopen("results.yaml", "w");
         fprintf(fp, "code: miniAMR\n");
         fprintf(fp, "version: %s\n", version);
         fprintf(fp, "ranks: %d\n", num_pes);
         fprintf(fp, "npx: %d\n", npx);
         fprintf(fp, "npy: %d\n", npy);
         fprintf(fp, "npz: %d\n", npz);
         fprintf(fp, "init_block_x: %d\n", init_block_x);
         fprintf(fp, "init_block_y: %d\n", init_block_y);
         fprintf(fp, "init_block_z: %d\n", init_block_z);
         fprintf(fp, "x_block_size: %d\n", x_block_size);
         fprintf(fp, "y_block_size: %d\n", y_block_size);
         fprintf(fp, "z_block_size: %d\n", z_block_size);
         fprintf(fp, "reorder: %d\n", reorder);
         fprintf(fp, "permute: %d\n", permute);
         fprintf(fp, "max_blocks_allowed: %d\n", max_num_blocks);
         fprintf(fp, "code: %d\n", code);
         fprintf(fp, "num_refine: %d\n", num_refine);
         fprintf(fp, "block_change: %d\n", block_change);
         fprintf(fp, "refine_ghost: %d\n", refine_ghost);
         fprintf(fp, "uniform_refine: %d\n", uniform_refine);
         fprintf(fp, "num_objects: %d\n", num_objects);
         for (i = 0; i < num_objects; i++) {
            op = &objects[i];
            fprintf(fp, "obj%dtype: %d\n", i, op->type);
            fprintf(fp, "obj%dbounce: %d\n", i, op->bounce);
            fprintf(fp, "obj%dcenter_x: %lf\n", i, op->cen[0]);
            fprintf(fp, "obj%dcenter_y: %lf\n", i, op->cen[1]);
            fprintf(fp, "obj%dcenter_z: %lf\n", i, op->cen[2]);
            fprintf(fp, "obj%dmove_x: %lf\n", i, op->move[0]);
            fprintf(fp, "obj%dmove_y: %lf\n", i, op->move[1]);
            fprintf(fp, "obj%dmove_z: %lf\n", i, op->move[2]);
            fprintf(fp, "obj%dsize_x: %lf\n", i, op->size[0]);
            fprintf(fp, "obj%dsize_y: %lf\n", i, op->size[1]);
            fprintf(fp, "obj%dsize_z: %lf\n", i, op->size[2]);
            fprintf(fp, "obj%dinc_x: %lf\n", i, op->inc[0]);
            fprintf(fp, "obj%dinc_y: %lf\n", i, op->inc[1]);
            fprintf(fp, "obj%dinc_z: %lf\n", i, op->inc[2]);
         }
         fprintf(fp, "num_tsteps: %d\n", num_tsteps);
         fprintf(fp, "stages_per_timestep: %d\n", stages_per_ts);
         fprintf(fp, "nonblocking: %d\n", nonblocking);
         fprintf(fp, "checksum_freq: %d\n", checksum_freq);
         fprintf(fp, "refine_freq: %d\n", refine_freq);
         fprintf(fp, "lb_opt: %d\n", lb_opt);
         fprintf(fp, "inbalance: %d\n", inbalance);
         fprintf(fp, "lb_method: %d\n", lb_method);
         fprintf(fp, "plot_freq: %d\n", plot_freq);
         fprintf(fp, "num_vars: %d\n", num_vars);
         fprintf(fp, "stencil: %d\n", stencil);
         fprintf(fp, "comm_vars: %d\n", comm_vars);
         fprintf(fp, "error_tol: %d\n", error_tol);
         fprintf(fp, "send_faces: %d\n", send_faces);

         fprintf(fp, "total_time_ave: %lf\n", average[0]);
         fprintf(fp, "total_time_min: %lf\n", minimum[0]);
         fprintf(fp, "total_time_max: %lf\n", maximum[0]);
         fprintf(fp, "memory_used_ave: %lf\n", average[112]);
         fprintf(fp, "memory_used_min: %lf\n", minimum[112]);
         fprintf(fp, "memory_used_max: %lf\n", maximum[112]);
         fprintf(fp, "compute_time_ave: %lf\n", average[38]);
         fprintf(fp, "compute_time_min: %lf\n", minimum[38]);
         fprintf(fp, "compute_time_max: %lf\n", maximum[38]);
         fprintf(fp, "total_gflops: %lf\n", total_gflops);
         fprintf(fp, "ave_gflops: %lf\n", gflops_rank);

         fprintf(fp, "total_comm_ave: %lf\n", average[37]);
         fprintf(fp, "total_comm_min: %lf\n", minimum[37]);
         fprintf(fp, "total_comm_max: %lf\n", maximum[37]);
         fprintf(fp, "   total_post_recv_ave: %lf\n", average[2]);
         fprintf(fp, "   total_post_recv_min: %lf\n", minimum[2]);
         fprintf(fp, "   total_post_recv_max: %lf\n", maximum[2]);
         fprintf(fp, "   total_pack_faces_ave: %lf\n", average[3]);
         fprintf(fp, "   total_pack_faces_min: %lf\n", minimum[3]);
         fprintf(fp, "   total_pack_faces_max: %lf\n", maximum[3]);
         fprintf(fp, "   total_send_mess_ave: %lf\n", average[4]);
         fprintf(fp, "   total_send_mess_min: %lf\n", minimum[4]);
         fprintf(fp, "   total_send_mess_max: %lf\n", maximum[4]);
         fprintf(fp, "   total_exch_same_ave: %lf\n", average[5]);
         fprintf(fp, "   total_exch_same_min: %lf\n", minimum[5]);
         fprintf(fp, "   total_exch_same_max: %lf\n", maximum[5]);
         fprintf(fp, "   total_exch_diff_ave: %lf\n", average[6]);
         fprintf(fp, "   total_exch_diff_min: %lf\n", minimum[6]);
         fprintf(fp, "   total_exch_diff_max: %lf\n", maximum[6]);
         fprintf(fp, "   total_apply_bc_ave: %lf\n", average[7]);
         fprintf(fp, "   total_apply_bc_min: %lf\n", minimum[7]);
         fprintf(fp, "   total_apply_bc_max: %lf\n", maximum[7]);
         fprintf(fp, "   total_wait_time_ave: %lf\n", average[8]);
         fprintf(fp, "   total_wait_time_min: %lf\n", minimum[8]);
         fprintf(fp, "   total_wait_time_max: %lf\n", maximum[8]);
         fprintf(fp, "   total_unpack_faces_ave: %lf\n", average[9]);
         fprintf(fp, "   total_unpack_faces_min: %lf\n", minimum[9]);
         fprintf(fp, "   total_unpack_faces_max: %lf\n", maximum[9]);
         fprintf(fp, "   comm_partners_total: %lf\n", average[120]);
         fprintf(fp, "   comm_partners_total_min: %lf\n", minimum[125]);
         fprintf(fp, "   comm_partners_total_max: %lf\n", maximum[130]);
         fprintf(fp, "   comm_partners_unique: %lf\n", average[121]);
         fprintf(fp, "   comm_partners_unique_min: %lf\n", minimum[126]);
         fprintf(fp, "   comm_partners_unique_max: %lf\n", maximum[131]);
         fprintf(fp, "   total_mess_recv_ave: %lf\n", average[71]);
         fprintf(fp, "   total_mess_recv_min: %lf\n", minimum[71]);
         fprintf(fp, "   total_mess_recv_max: %lf\n", maximum[71]);
         fprintf(fp, "   total_byte_recv_ave: %lf\n", average[69]);
         fprintf(fp, "   total_byte_recv_min: %lf\n", minimum[69]);
         fprintf(fp, "   total_byte_recv_max: %lf\n", maximum[69]);
         fprintf(fp, "   total_face_recv_ave: %lf\n", average[73]);
         fprintf(fp, "   total_face_recv_min: %lf\n", minimum[73]);
         fprintf(fp, "   total_face_recv_max: %lf\n", maximum[73]);
         fprintf(fp, "   total_mess_send_ave: %lf\n", average[72]);
         fprintf(fp, "   total_mess_send_min: %lf\n", minimum[72]);
         fprintf(fp, "   total_mess_send_max: %lf\n", maximum[72]);
         fprintf(fp, "   total_byte_send_ave: %lf\n", average[70]);
         fprintf(fp, "   total_byte_send_min: %lf\n", minimum[70]);
         fprintf(fp, "   total_byte_send_max: %lf\n", maximum[70]);
         fprintf(fp, "   total_face_send_ave: %lf\n", average[74]);
         fprintf(fp, "   total_face_send_min: %lf\n", minimum[74]);
         fprintf(fp, "   total_face_send_max: %lf\n", maximum[74]);
         fprintf(fp, "   total_face_exch_same_ave: %lf\n", average[76]);
         fprintf(fp, "   total_face_exch_same_min: %lf\n", minimum[76]);
         fprintf(fp, "   total_face_exch_same_max: %lf\n", maximum[76]);
         fprintf(fp, "   total_face_exch_diff_ave: %lf\n", average[77]);
         fprintf(fp, "   total_face_exch_diff_min: %lf\n", minimum[77]);
         fprintf(fp, "   total_face_exch_diff_max: %lf\n", maximum[77]);
         fprintf(fp, "   total_face_bc_apply_ave: %lf\n", average[75]);
         fprintf(fp, "   total_face_bc_apply_min: %lf\n", minimum[75]);
         fprintf(fp, "   total_face_bc_apply_max: %lf\n", maximum[75]);

         fprintf(fp, "   x_comm_ave: %lf\n", average[10]);
         fprintf(fp, "   x_comm_min: %lf\n", minimum[10]);
         fprintf(fp, "   x_comm_max: %lf\n", maximum[10]);
         fprintf(fp, "      x_post_recv_ave: %lf\n", average[11]);
         fprintf(fp, "      x_post_recv_min: %lf\n", minimum[11]);
         fprintf(fp, "      x_post_recv_max: %lf\n", maximum[11]);
         fprintf(fp, "      x_pack_faces_ave: %lf\n", average[12]);
         fprintf(fp, "      x_pack_faces_min: %lf\n", minimum[12]);
         fprintf(fp, "      x_pack_faces_max: %lf\n", maximum[12]);
         fprintf(fp, "      x_send_mess_ave: %lf\n", average[13]);
         fprintf(fp, "      x_send_mess_min: %lf\n", minimum[13]);
         fprintf(fp, "      x_send_mess_max: %lf\n", maximum[13]);
         fprintf(fp, "      x_exch_same_ave: %lf\n", average[14]);
         fprintf(fp, "      x_exch_same_min: %lf\n", minimum[14]);
         fprintf(fp, "      x_exch_same_max: %lf\n", maximum[14]);
         fprintf(fp, "      x_exch_diff_ave: %lf\n", average[15]);
         fprintf(fp, "      x_exch_diff_min: %lf\n", minimum[15]);
         fprintf(fp, "      x_exch_diff_max: %lf\n", maximum[15]);
         fprintf(fp, "      x_apply_bc_ave: %lf\n", average[16]);
         fprintf(fp, "      x_apply_bc_min: %lf\n", minimum[16]);
         fprintf(fp, "      x_apply_bc_max: %lf\n", maximum[16]);
         fprintf(fp, "      x_wait_time_ave: %lf\n", average[17]);
         fprintf(fp, "      x_wait_time_min: %lf\n", minimum[17]);
         fprintf(fp, "      x_wait_time_max: %lf\n", maximum[17]);
         fprintf(fp, "      x_unpack_faces_ave: %lf\n", average[18]);
         fprintf(fp, "      x_unpack_faces_min: %lf\n", minimum[18]);
         fprintf(fp, "      x_unpack_faces_max: %lf\n", maximum[18]);
         fprintf(fp, "      x_comm_partners: %lf\n", average[117]);
         fprintf(fp, "      x_comm_partners_min: %d\n", (int) minimum[122]);
         fprintf(fp, "      x_comm_partners_max: %d\n", (int) maximum[127]);
         fprintf(fp, "      x_mess_recv_ave: %lf\n", average[80]);
         fprintf(fp, "      x_mess_recv_min: %lf\n", minimum[80]);
         fprintf(fp, "      x_mess_recv_max: %lf\n", maximum[80]);
         fprintf(fp, "      x_byte_recv_ave: %lf\n", average[78]);
         fprintf(fp, "      x_byte_recv_min: %lf\n", minimum[78]);
         fprintf(fp, "      x_byte_recv_max: %lf\n", maximum[78]);
         fprintf(fp, "      x_face_recv_ave: %lf\n", average[82]);
         fprintf(fp, "      x_face_recv_min: %lf\n", minimum[82]);
         fprintf(fp, "      x_face_recv_max: %lf\n", maximum[82]);
         fprintf(fp, "      x_mess_send_ave: %lf\n", average[81]);
         fprintf(fp, "      x_mess_send_min: %lf\n", minimum[81]);
         fprintf(fp, "      x_mess_send_max: %lf\n", maximum[81]);
         fprintf(fp, "      x_byte_send_ave: %lf\n", average[79]);
         fprintf(fp, "      x_byte_send_min: %lf\n", minimum[79]);
         fprintf(fp, "      x_byte_send_max: %lf\n", maximum[79]);
         fprintf(fp, "      x_face_send_ave: %lf\n", average[83]);
         fprintf(fp, "      x_face_send_min: %lf\n", minimum[83]);
         fprintf(fp, "      x_face_send_max: %lf\n", maximum[83]);
         fprintf(fp, "      x_face_exch_same_ave: %lf\n", average[85]);
         fprintf(fp, "      x_face_exch_same_min: %lf\n", minimum[85]);
         fprintf(fp, "      x_face_exch_same_max: %lf\n", maximum[85]);
         fprintf(fp, "      x_face_exch_diff_ave: %lf\n", average[86]);
         fprintf(fp, "      x_face_exch_diff_min: %lf\n", minimum[86]);
         fprintf(fp, "      x_face_exch_diff_max: %lf\n", maximum[86]);
         fprintf(fp, "      x_face_bc_apply_ave: %lf\n", average[84]);
         fprintf(fp, "      x_face_bc_apply_min: %lf\n", minimum[84]);
         fprintf(fp, "      x_face_bc_apply_max: %lf\n", maximum[84]);

         fprintf(fp, "   y_comm_ave: %lf\n", average[19]);
         fprintf(fp, "   y_comm_min: %lf\n", minimum[19]);
         fprintf(fp, "   y_comm_max: %lf\n", maximum[19]);
         fprintf(fp, "      y_post_recv_ave: %lf\n", average[20]);
         fprintf(fp, "      y_post_recv_min: %lf\n", minimum[20]);
         fprintf(fp, "      y_post_recv_max: %lf\n", maximum[20]);
         fprintf(fp, "      y_pack_faces_ave: %lf\n", average[21]);
         fprintf(fp, "      y_pack_faces_min: %lf\n", minimum[21]);
         fprintf(fp, "      y_pack_faces_max: %lf\n", maximum[21]);
         fprintf(fp, "      y_send_mess_ave: %lf\n", average[22]);
         fprintf(fp, "      y_send_mess_min: %lf\n", minimum[22]);
         fprintf(fp, "      y_send_mess_max: %lf\n", maximum[22]);
         fprintf(fp, "      y_exch_same_ave: %lf\n", average[23]);
         fprintf(fp, "      y_exch_same_min: %lf\n", minimum[23]);
         fprintf(fp, "      y_exch_same_max: %lf\n", maximum[23]);
         fprintf(fp, "      y_exch_diff_ave: %lf\n", average[24]);
         fprintf(fp, "      y_exch_diff_min: %lf\n", minimum[24]);
         fprintf(fp, "      y_exch_diff_max: %lf\n", maximum[24]);
         fprintf(fp, "      y_apply_bc_ave: %lf\n", average[25]);
         fprintf(fp, "      y_apply_bc_min: %lf\n", minimum[25]);
         fprintf(fp, "      y_apply_bc_max: %lf\n", maximum[25]);
         fprintf(fp, "      y_wait_time_ave: %lf\n", average[26]);
         fprintf(fp, "      y_wait_time_min: %lf\n", minimum[26]);
         fprintf(fp, "      y_wait_time_max: %lf\n", maximum[26]);
         fprintf(fp, "      y_unpack_faces_ave: %lf\n", average[27]);
         fprintf(fp, "      y_unpack_faces_min: %lf\n", minimum[27]);
         fprintf(fp, "      y_unpack_faces_max: %lf\n", maximum[27]);
         fprintf(fp, "      y_comm_partners: %lf\n", average[118]);
         fprintf(fp, "      y_comm_partners_min: %d\n", (int) minimum[123]);
         fprintf(fp, "      y_comm_partners_max: %d\n", (int) maximum[128]);
         fprintf(fp, "      y_mess_recv_ave: %lf\n", average[89]);
         fprintf(fp, "      y_mess_recv_min: %lf\n", minimum[89]);
         fprintf(fp, "      y_mess_recv_max: %lf\n", maximum[89]);
         fprintf(fp, "      y_byte_recv_ave: %lf\n", average[87]);
         fprintf(fp, "      y_byte_recv_min: %lf\n", minimum[87]);
         fprintf(fp, "      y_byte_recv_max: %lf\n", maximum[87]);
         fprintf(fp, "      y_face_recv_ave: %lf\n", average[91]);
         fprintf(fp, "      y_face_recv_min: %lf\n", minimum[91]);
         fprintf(fp, "      y_face_recv_max: %lf\n", maximum[91]);
         fprintf(fp, "      y_mess_send_ave: %lf\n", average[90]);
         fprintf(fp, "      y_mess_send_min: %lf\n", minimum[90]);
         fprintf(fp, "      y_mess_send_max: %lf\n", maximum[90]);
         fprintf(fp, "      y_byte_send_ave: %lf\n", average[88]);
         fprintf(fp, "      y_byte_send_min: %lf\n", minimum[88]);
         fprintf(fp, "      y_byte_send_max: %lf\n", maximum[88]);
         fprintf(fp, "      y_face_send_ave: %lf\n", average[92]);
         fprintf(fp, "      y_face_send_min: %lf\n", minimum[92]);
         fprintf(fp, "      y_face_send_max: %lf\n", maximum[92]);
         fprintf(fp, "      y_face_exch_same_ave: %lf\n", average[94]);
         fprintf(fp, "      y_face_exch_same_min: %lf\n", minimum[94]);
         fprintf(fp, "      y_face_exch_same_max: %lf\n", maximum[94]);
         fprintf(fp, "      y_face_exch_diff_ave: %lf\n", average[95]);
         fprintf(fp, "      y_face_exch_diff_min: %lf\n", minimum[95]);
         fprintf(fp, "      y_face_exch_diff_max: %lf\n", maximum[95]);
         fprintf(fp, "      y_face_bc_apply_ave: %lf\n", average[93]);
         fprintf(fp, "      y_face_bc_apply_min: %lf\n", minimum[93]);
         fprintf(fp, "      y_face_bc_apply_max: %lf\n", maximum[93]);

         fprintf(fp, "   z_comm_ave: %lf\n", average[28]);
         fprintf(fp, "   z_comm_min: %lf\n", minimum[28]);
         fprintf(fp, "   z_comm_max: %lf\n", maximum[28]);
         fprintf(fp, "      z_post_recv_ave: %lf\n", average[29]);
         fprintf(fp, "      z_post_recv_min: %lf\n", minimum[29]);
         fprintf(fp, "      z_post_recv_max: %lf\n", maximum[29]);
         fprintf(fp, "      z_pack_faces_ave: %lf\n", average[30]);
         fprintf(fp, "      z_pack_faces_min: %lf\n", minimum[30]);
         fprintf(fp, "      z_pack_faces_max: %lf\n", maximum[30]);
         fprintf(fp, "      z_send_mess_ave: %lf\n", average[31]);
         fprintf(fp, "      z_send_mess_min: %lf\n", minimum[31]);
         fprintf(fp, "      z_send_mess_max: %lf\n", maximum[31]);
         fprintf(fp, "      z_exch_same_ave: %lf\n", average[32]);
         fprintf(fp, "      z_exch_same_min: %lf\n", minimum[32]);
         fprintf(fp, "      z_exch_same_max: %lf\n", maximum[32]);
         fprintf(fp, "      z_exch_diff_ave: %lf\n", average[33]);
         fprintf(fp, "      z_exch_diff_min: %lf\n", minimum[33]);
         fprintf(fp, "      z_exch_diff_max: %lf\n", maximum[33]);
         fprintf(fp, "      z_apply_bc_ave: %lf\n", average[34]);
         fprintf(fp, "      z_apply_bc_min: %lf\n", minimum[34]);
         fprintf(fp, "      z_apply_bc_max: %lf\n", maximum[34]);
         fprintf(fp, "      z_wait_time_ave: %lf\n", average[35]);
         fprintf(fp, "      z_wait_time_min: %lf\n", minimum[35]);
         fprintf(fp, "      z_wait_time_max: %lf\n", maximum[35]);
         fprintf(fp, "      z_unpack_faces_ave: %lf\n", average[36]);
         fprintf(fp, "      z_unpack_faces_min: %lf\n", minimum[36]);
         fprintf(fp, "      z_unpack_faces_max: %lf\n", maximum[36]);
         fprintf(fp, "      z_comm_partners: %lf\n", average[119]);
         fprintf(fp, "      z_comm_partners_min: %d\n", (int) minimum[124]);
         fprintf(fp, "      z_comm_partners_max: %d\n", (int) maximum[129]);
         fprintf(fp, "      z_mess_recv_ave: %lf\n", average[98]);
         fprintf(fp, "      z_mess_recv_min: %lf\n", minimum[98]);
         fprintf(fp, "      z_mess_recv_max: %lf\n", maximum[98]);
         fprintf(fp, "      z_byte_recv_ave: %lf\n", average[96]);
         fprintf(fp, "      z_byte_recv_min: %lf\n", minimum[96]);
         fprintf(fp, "      z_byte_recv_max: %lf\n", maximum[96]);
         fprintf(fp, "      z_face_recv_ave: %lf\n", average[100]);
         fprintf(fp, "      z_face_recv_min: %lf\n", minimum[100]);
         fprintf(fp, "      z_face_recv_max: %lf\n", maximum[100]);
         fprintf(fp, "      z_mess_send_ave: %lf\n", average[99]);
         fprintf(fp, "      z_mess_send_min: %lf\n", minimum[99]);
         fprintf(fp, "      z_mess_send_max: %lf\n", maximum[99]);
         fprintf(fp, "      z_byte_send_ave: %lf\n", average[97]);
         fprintf(fp, "      z_byte_send_min: %lf\n", minimum[97]);
         fprintf(fp, "      z_byte_send_max: %lf\n", maximum[97]);
         fprintf(fp, "      z_face_send_ave: %lf\n", average[101]);
         fprintf(fp, "      z_face_send_min: %lf\n", minimum[101]);
         fprintf(fp, "      z_face_send_max: %lf\n", maximum[101]);
         fprintf(fp, "      z_face_exch_same_ave: %lf\n", average[103]);
         fprintf(fp, "      z_face_exch_same_min: %lf\n", minimum[103]);
         fprintf(fp, "      z_face_exch_same_max: %lf\n", maximum[103]);
         fprintf(fp, "      z_face_exch_diff_ave: %lf\n", average[104]);
         fprintf(fp, "      z_face_exch_diff_min: %lf\n", minimum[104]);
         fprintf(fp, "      z_face_exch_diff_max: %lf\n", maximum[104]);
         fprintf(fp, "      z_face_bc_apply_ave: %lf\n", average[102]);
         fprintf(fp, "      z_face_bc_apply_min: %lf\n", minimum[102]);
         fprintf(fp, "      z_face_bc_apply_max: %lf\n", maximum[102]);

         fprintf(fp, "gridsum_time_ave: %lf\n", average[39]);
         fprintf(fp, "gridsum_time_min: %lf\n", minimum[39]);
         fprintf(fp, "gridsum_time_max: %lf\n", maximum[39]);
         fprintf(fp, "   gridsum_reduce_ave: %lf\n", average[40]);
         fprintf(fp, "   gridsum_reduce_min: %lf\n", minimum[40]);
         fprintf(fp, "   gridsum_reduce_max: %lf\n", maximum[40]);
         fprintf(fp, "   gridsum_calc_ave: %lf\n", average[41]);
         fprintf(fp, "   gridsum_calc_min: %lf\n", minimum[41]);
         fprintf(fp, "   gridsum_calc_max: %lf\n", maximum[41]);

         fprintf(fp, "refine_time_ave: %lf\n", average[42]);
         fprintf(fp, "refine_time_min: %lf\n", minimum[42]);
         fprintf(fp, "refine_time_max: %lf\n", maximum[42]);
         fprintf(fp, "   total_blocks_ts_ave: %lf\n",
                ((double) total_blocks)/((double) (num_tsteps*stages_per_ts)));
         fprintf(fp, "   total_blocks_ts_min: %lld\n", (long long) nb_min);
         fprintf(fp, "   total_blocks_ts_max: %lld\n", (long long) nb_max);
         fprintf(fp, "   blocks_split_ave: %lf\n", average[105]);
         fprintf(fp, "   blocks_split_min: %lf\n", minimum[105]);
         fprintf(fp, "   blocks_split_max: %lf\n", maximum[105]);
         fprintf(fp, "   blocks_reformed_ave: %lf\n", average[106]);
         fprintf(fp, "   blocks_reformed_min: %lf\n", minimum[106]);
         fprintf(fp, "   blocks_reformed_max: %lf\n", maximum[106]);
         fprintf(fp, "   blocks_moved_tot_ave: %lf\n", average[107]);
         fprintf(fp, "   blocks_moved_tot_min: %lf\n", minimum[107]);
         fprintf(fp, "   blocks_moved_tot_max: %lf\n", maximum[107]);
         fprintf(fp, "   blocks_moved_lb_ave: %lf\n", average[108]);
         fprintf(fp, "   blocks_moved_lb_min: %lf\n", minimum[108]);
         fprintf(fp, "   blocks_moved_lb_max: %lf\n", maximum[108]);
         fprintf(fp, "   blocks_moved_redist_ave: %lf\n", average[109]);
         fprintf(fp, "   blocks_moved_redist_min: %lf\n", minimum[109]);
         fprintf(fp, "   blocks_moved_redist_max: %lf\n", maximum[109]);
         fprintf(fp, "   blocks_moved_coarsen_ave: %lf\n", average[110]);
         fprintf(fp, "   blocks_moved_coarsen_min: %lf\n", minimum[110]);
         fprintf(fp, "   blocks_moved_coarsen_max: %lf\n", maximum[110]);
         fprintf(fp, "   time_compare_obj_ave: %lf\n", average[43]);
         fprintf(fp, "   time_compare_obj_min: %lf\n", minimum[43]);
         fprintf(fp, "   time_compare_obj_max: %lf\n", maximum[43]);
         fprintf(fp, "   time_mark_refine_ave: %lf\n", average[44]);
         fprintf(fp, "   time_mark_refine_min: %lf\n", minimum[44]);
         fprintf(fp, "   time_mark_refine_max: %lf\n", maximum[44]);
         fprintf(fp, "   time_comm_block1_ave: %lf\n", average[47]);
         fprintf(fp, "   time_comm_block1_min: %lf\n", minimum[47]);
         fprintf(fp, "   time_comm_block1_max: %lf\n", maximum[47]);
         fprintf(fp, "   time_split_block_ave: %lf\n", average[46]);
         fprintf(fp, "   time_split_block_min: %lf\n", minimum[46]);
         fprintf(fp, "   time_split_block_max: %lf\n", maximum[46]);
         fprintf(fp, "   time_comm_block2_ave: %lf\n", average[48]);
         fprintf(fp, "   time_comm_block2_min: %lf\n", minimum[48]);
         fprintf(fp, "   time_comm_block2_max: %lf\n", maximum[48]);
         fprintf(fp, "   time_sync_ave: %lf\n", average[49]);
         fprintf(fp, "   time_sync_min: %lf\n", minimum[49]);
         fprintf(fp, "   time_sync_max: %lf\n", maximum[49]);
         fprintf(fp, "   time_misc_ave: %lf\n", average[45]);
         fprintf(fp, "   time_misc_min: %lf\n", minimum[45]);
         fprintf(fp, "   time_misc_max: %lf\n", maximum[45]);
         fprintf(fp, "   time_total_coarsen_ave: %lf\n", average[50]);
         fprintf(fp, "   time_total_coarsen_min: %lf\n", minimum[50]);
         fprintf(fp, "   time_total_coarsen_max: %lf\n", maximum[50]);
         fprintf(fp, "      time_coarsen_ave: %lf\n", average[51]);
         fprintf(fp, "      time_coarsen_min: %lf\n", minimum[51]);
         fprintf(fp, "      time_coarsen_max: %lf\n", maximum[51]);
         fprintf(fp, "      time_coarsen_pack_ave: %lf\n", average[52]);
         fprintf(fp, "      time_coarsen_pack_min: %lf\n", minimum[52]);
         fprintf(fp, "      time_coarsen_pack_max: %lf\n", maximum[52]);
         fprintf(fp, "      time_coarsen_move_ave: %lf\n", average[53]);
         fprintf(fp, "      time_coarsen_move_min: %lf\n", minimum[53]);
         fprintf(fp, "      time_coarsen_move_max: %lf\n", maximum[53]);
         fprintf(fp, "      time_coarsen_unpack_ave: %lf\n", average[54]);
         fprintf(fp, "      time_coarsen_unpack_min: %lf\n", minimum[54]);
         fprintf(fp, "      time_coarsen_unpack_max: %lf\n", maximum[54]);
         fprintf(fp, "   time_total_redist_ave: %lf\n", average[64]);
         fprintf(fp, "   time_total_redist_min: %lf\n", minimum[64]);
         fprintf(fp, "   time_total_redist_max: %lf\n", maximum[64]);
         fprintf(fp, "      time_redist_choose_ave: %lf\n", average[65]);
         fprintf(fp, "      time_redist_choose_min: %lf\n", minimum[65]);
         fprintf(fp, "      time_redist_choose_max: %lf\n", maximum[65]);
         fprintf(fp, "      time_redist_pack_ave: %lf\n", average[66]);
         fprintf(fp, "      time_redist_pack_min: %lf\n", minimum[66]);
         fprintf(fp, "      time_redist_pack_max: %lf\n", maximum[66]);
         fprintf(fp, "      time_redist_move_ave: %lf\n", average[67]);
         fprintf(fp, "      time_redist_move_min: %lf\n", minimum[67]);
         fprintf(fp, "      time_redist_move_max: %lf\n", maximum[67]);
         fprintf(fp, "      time_redist_unpack_ave: %lf\n", average[68]);
         fprintf(fp, "      time_redist_unpack_min: %lf\n", minimum[68]);
         fprintf(fp, "      time_redist_unpack_max: %lf\n", maximum[68]);
         fprintf(fp, "   time_total_load_bal_ave: %lf\n", average[55]);
         fprintf(fp, "   time_total_load_bal_min: %lf\n", minimum[55]);
         fprintf(fp, "   time_total_load_bal_max: %lf\n", maximum[55]);
         fprintf(fp, "      time_load_bal_sort_ave: %lf\n", average[56]);
         fprintf(fp, "      time_load_bal_sort_min: %lf\n", minimum[56]);
         fprintf(fp, "      time_load_bal_sort_max: %lf\n", maximum[56]);
         fprintf(fp, "      time_lb_move_dots_ave: %lf\n", average[61]);
         fprintf(fp, "      time_lb_move_dots_min: %lf\n", minimum[61]);
         fprintf(fp, "      time_lb_move_dots_max: %lf\n", maximum[61]);
         fprintf(fp, "      time_lb_move_blocks_ave: %lf\n", average[62]);
         fprintf(fp, "      time_lb_move_blocks_min: %lf\n", minimum[62]);
         fprintf(fp, "      time_lb_move_blocks_max: %lf\n", maximum[62]);
         fprintf(fp, "         time_lb_mb_pack_ave: %lf\n", average[57]);
         fprintf(fp, "         time_lb_mb_pack_min: %lf\n", minimum[57]);
         fprintf(fp, "         time_lb_mb_pack_max: %lf\n", maximum[57]);
         fprintf(fp, "         time_lb_mb_move_ave: %lf\n", average[58]);
         fprintf(fp, "         time_lb_mb_move_min: %lf\n", minimum[58]);
         fprintf(fp, "         time_lb_mb_move_max: %lf\n", maximum[58]);
         fprintf(fp, "         time_lb_mb_unpack_ave: %lf\n", average[59]);
         fprintf(fp, "         time_lb_mb_unpack_min: %lf\n", minimum[59]);
         fprintf(fp, "         time_lb_mb_unpack_max: %lf\n", maximum[59]);
         fprintf(fp, "         time_lb_mb_misc_ave: %lf\n", average[60]);
         fprintf(fp, "         time_lb_mb_misc_min: %lf\n", minimum[60]);
         fprintf(fp, "         time_lb_mb_misc_max: %lf\n", maximum[60]);

         fprintf(fp, "plot_time_ave: %lf\n", average[63]);
         fprintf(fp, "plot_time_min: %lf\n", minimum[63]);
         fprintf(fp, "plot_time_max: %lf\n", maximum[63]);

         fclose(fp);
      }

      if (report_perf & 2) {
         fp = fopen("results.txt", "w");

         fprintf(fp, "\n ================ Start report ===================\n\n");
         fprintf(fp, "          Mantevo miniAMR\n");
         fprintf(fp, "          version %s\n\n", version);

         fprintf(fp, "Run on %d ranks arranged in a %d x %d x %d grid\n", num_pes,
                npx, npy, npz);
         fprintf(fp, "initial blocks per rank %d x %d x %d\n", init_block_x,
                init_block_y, init_block_z);
         fprintf(fp, "block size %d x %d x %d\n", x_block_size, y_block_size,
                z_block_size);
         if (reorder)
            fprintf(fp, "Initial ranks arranged by RCB across machine\n\n");
         else
            fprintf(fp, "Initial ranks arranged as a grid across machine\n\n");
         if (permute)
            fprintf(fp, "Order of exchanges permuted\n");
         fprintf(fp, "Maximum number of blocks per rank is %d\n",
                 max_num_blocks);
         if (code)
            fprintf(fp, "Code set to code %d\n", code);
         fprintf(fp, "Number of levels of refinement is %d\n", num_refine);
         fprintf(fp, "Blocks can change by %d levels per refinement step\n",
            block_change);
         if (refine_ghost)
            fprintf(fp, "Ghost cells will be used determine is block is refined\n");
         if (uniform_refine)
            fprintf(fp, "\nBlocks will be uniformly refined\n");
         else {
            fprintf(fp, "\nBlocks will be refined by %d objects\n\n", num_objects);
            for (i = 0; i < num_objects; i++) {
               op = &objects[i];
               if (op->type == 0)
                  fprintf(fp, "Object %d is the surface of a rectangle\n", i);
               else if (op->type == 1)
                  fprintf(fp, "Object %d is the volume of a rectangle\n", i);
               else if (op->type == 2)
                  fprintf(fp, "Object %d is the surface of a spheroid\n", i);
               else if (op->type == 3)
                  fprintf(fp, "Object %d is the volume of a spheroid\n", i);
               else if (op->type == 4)
                  fprintf(fp, "Object %d is the surface of x+ hemispheroid\n", i);
               else if (op->type == 5)
                  fprintf(fp, "Object %d is the volume of x+ hemispheroid\n", i);
               else if (op->type == 6)
                  fprintf(fp, "Object %d is the surface of x- hemispheroid\n", i);
               else if (op->type == 7)
                  fprintf(fp, "Object %d is the volume of x- hemispheroid\n", i);
               else if (op->type == 8)
                  fprintf(fp, "Object %d is the surface of y+ hemispheroid\n", i);
               else if (op->type == 9)
                  fprintf(fp, "Object %d is the volume of y+ hemispheroid\n", i);
               else if (op->type == 10)
                  fprintf(fp, "Object %d is the surface of y- hemispheroid\n", i);
               else if (op->type == 11)
                  fprintf(fp, "Object %d is the volume of y- hemispheroid\n", i);
               else if (op->type == 12)
                  fprintf(fp, "Object %d is the surface of z+ hemispheroid\n", i);
               else if (op->type == 13)
                  fprintf(fp, "Object %d is the volume of z+ hemispheroid\n", i);
               else if (op->type == 14)
                  fprintf(fp, "Object %d is the surface of z- hemispheroid\n", i);
               else if (op->type == 15)
                  fprintf(fp, "Object %d is the volume of z- hemispheroid\n", i);
               else if (op->type == 20)
                  fprintf(fp, "Object %d is the surface of x axis cylinder\n", i);
               else if (op->type == 21)
                  fprintf(fp, "Object %d is the volune of x axis cylinder\n", i);
               else if (op->type == 22)
                  fprintf(fp, "Object %d is the surface of y axis cylinder\n", i);
               else if (op->type == 23)
                  fprintf(fp, "Object %d is the volune of y axis cylinder\n", i);
               else if (op->type == 24)
                  fprintf(fp, "Object %d is the surface of z axis cylinder\n", i);
               else if (op->type == 25)
                  fprintf(fp, "Object %d is the volune of z axis cylinder\n", i);
               if (op->bounce == 0)
                  fprintf(fp, "Object may leave mesh\n");
               else
                  fprintf(fp, "Object center will bounce off of walls\n");
               fprintf(fp, "Center starting at %lf %lf %lf\n",
                      op->orig_cen[0], op->orig_cen[1], op->orig_cen[2]);
               fprintf(fp, "Center end at %lf %lf %lf\n",
                      op->cen[0], op->cen[1], op->cen[2]);

               if (use_time) {
                  fprintf(fp, "Velocity of %lf %lf %lf\n",
                         op->orig_move[0], op->orig_move[1], op->orig_move[2]);
                  delta = end_time/((double) num_tsteps);
                  fprintf(fp,
                         "   Rate relative to smallest cell size %lf %lf %lf\n",
                         delta*op->orig_move[0]*((double)
                                                 (mesh_size[0]*x_block_size)),
                         delta*op->orig_move[1]*((double)
                                                 (mesh_size[1]*y_block_size)),
                         delta*op->orig_move[2]*((double)
                                                 (mesh_size[2]*z_block_size)));
               } else {
                  fprintf(fp, "Moving at %lf %lf %lf per timestep\n",
                         op->orig_move[0], op->orig_move[1], op->orig_move[2]);
                  fprintf(fp,
                      "   Rate relative to smallest cell size %lf %lf %lf\n",
                      op->orig_move[0]*((double) (mesh_size[0]*x_block_size)),
                      op->orig_move[1]*((double) (mesh_size[1]*y_block_size)),
                      op->orig_move[2]*((double) (mesh_size[2]*z_block_size)));
               }
               fprintf(fp, "Initial size %lf %lf %lf\n",
                      op->orig_size[0], op->orig_size[1], op->orig_size[2]);
               fprintf(fp, "Final size %lf %lf %lf\n",
                      op->size[0], op->size[1], op->size[2]);
               fprintf(fp, "Size increasing %lf %lf %lf per timestep\n",
                      op->inc[0], op->inc[1], op->inc[2]);
               fprintf(fp, "   Rate relative to smallest cell size %lf %lf %lf\n\n",
                      op->inc[0]*((double) (mesh_size[0]*x_block_size)),
                      op->inc[1]*((double) (mesh_size[1]*y_block_size)),
                      op->inc[2]*((double) (mesh_size[2]*z_block_size)));
            }
         }
         if (use_time)
            fprintf(fp, "\nTime %lf in %d timesteps\n", end_time, num_tsteps);
         else
            fprintf(fp, "\nNumber of timesteps is %d\n", num_tsteps);
         fprintf(fp, "Communicaion/computation stages per timestep is %d\n",
                stages_per_ts);
         if (nonblocking)
            fprintf(fp, "Communication will be performed with nonblocking sends\n");
         else
            fprintf(fp, "Communication will be performed with blocking sends\n");
         fprintf(fp, "Will perform checksums every %d stages\n", checksum_freq);
         fprintf(fp, "Will refine every %d timesteps\n", refine_freq);
         if (!lb_method)
            fprintf(fp, "Load balance by RCB (Recursive Coordinate Bisection)\n");
         else if (lb_method == 1)
            fprintf(fp, "Load balance by Morton Space Filling Curve\n");
         else if (lb_method == 2)
            fprintf(fp, "Load balance by Peano Hilbert Space Filling Curve\n");
         else
            fprintf(fp, "Load balance by truncated Hilbert Space Filling Curve\n");
         if (lb_opt == 0)
            fprintf(fp, "Load balance will not be performed\n");
         else
            fprintf(fp, "Load balance when inbalanced by %d%%\n", inbalance);
         if (lb_opt == 2)
            fprintf(fp, "Load balance at each phase of refinement step\n");
         if (plot_freq)
            fprintf(fp, "Will plot results every %d timesteps\n", plot_freq);
         else
            fprintf(fp, "Will not plot results\n");
         if (stencil)
            fprintf(fp, "Calculate on %d variables with %d point stencil\n",
                    num_vars, stencil);
         else
            fprintf(fp, "Calculate on %d variables with variable stencils\n",
                    num_vars);
         fprintf(fp, "Communicate %d variables at a time\n", comm_vars);
         fprintf(fp, "Error tolorance for variable sums is 10^(-%d)\n", error_tol);
         if (send_faces)
            fprintf(fp, "Will send data from each face seperately\n");

         fprintf(fp, "\nTotal time for test: ave, std, min, max (sec): %lf %lf %lf %lf\n\n",
                average[0], stddev[0], minimum[0], maximum[0]);

         fprintf(fp, "\nNumber of malloc calls: ave, std, min, max (sec): %lf %lf %lf %lf\n",
                average[111], stddev[111], minimum[111], maximum[111]);
         fprintf(fp, "\nAmount malloced: ave, std, min, max: %lf %lf %lf %lf\n",
                average[112], stddev[112], minimum[112], maximum[112]);
         fprintf(fp, "\nMalloc calls in init: ave, std, min, max (sec): %lf %lf %lf %lf\n",
                average[113], stddev[113], minimum[113], maximum[113]);
         fprintf(fp, "\nAmount malloced in init: ave, std, min, max: %lf %lf %lf %lf\n",
                average[114], stddev[114], minimum[114], maximum[114]);
         fprintf(fp, "\nMalloc calls in timestepping: ave, std, min, max (sec): %lf %lf %lf %lf\n",
                average[115], stddev[115], minimum[115], maximum[115]);
         fprintf(fp, "\nAmount malloced in timestepping: ave, std, min, max: %lf %lf %lf %lf\n\n",
                average[116], stddev[116], minimum[116], maximum[116]);

         fprintf(fp, "---------------------------------------------\n");
         fprintf(fp, "          Computational Performance\n");
         fprintf(fp, "---------------------------------------------\n\n");
         fprintf(fp, "     Time: ave, stddev, min, max (sec): %lf %lf %lf %lf\n\n",
                average[38], stddev[38], minimum[38], maximum[38]);
         fprintf(fp, "     total GFLOPS:             %lf\n", total_gflops);
         fprintf(fp, "     Average GFLOPS per rank:  %lf\n\n", gflops_rank);
         fprintf(fp, "     Total floating point ops: %lf\n\n", total_fp_ops);
         fprintf(fp, "        Adds:                  %lf\n", average[139]);
         fprintf(fp, "        Muls:                  %lf\n", average[140]);
         fprintf(fp, "        Divides:               %lf\n\n", average[141]);

         fprintf(fp, "---------------------------------------------\n");
         fprintf(fp, "           Interblock communication\n");
         fprintf(fp, "---------------------------------------------\n\n");
         fprintf(fp, "     Time: ave, stddev, min, max (sec): %lf %lf %lf %lf\n\n",
                average[37], stddev[37], minimum[37], maximum[37]);
         for (i = 0; i < 4; i++) {
            if (i == 0)
               fprintf(fp, "\nTotal communication:\n\n");
            else if (i == 1)
               fprintf(fp, "\nX direction communication statistics:\n\n");
            else if (i == 2)
               fprintf(fp, "\nY direction communication statistics:\n\n");
            else
               fprintf(fp, "\nZ direction communication statistics:\n\n");
            fprintf(fp, "                              average    stddev  minimum  maximum\n");
            fprintf(fp, "     Total                  : %lf %lf %lf %lf\n",
                   average[1+9*i], stddev[1+9*i], minimum[1+9*i],
                   maximum[1+9*i]);
            fprintf(fp, "     Post IRecv             : %lf %lf %lf %lf\n",
                   average[2+9*i], stddev[2+9*i], minimum[2+9*i],
                   maximum[2+9*i]);
            fprintf(fp, "     Pack faces             : %lf %lf %lf %lf\n",
                   average[3+9*i], stddev[3+9*i], minimum[3+9*i],
                   maximum[3+9*i]);
            fprintf(fp, "     Send messages          : %lf %lf %lf %lf\n",
                   average[4+9*i], stddev[4+9*i], minimum[4+9*i],
                   maximum[4+9*i]);
            fprintf(fp, "     Exchange same level    : %lf %lf %lf %lf\n",
                   average[5+9*i], stddev[5+9*i], minimum[5+9*i],
                   maximum[5+9*i]);
            fprintf(fp, "     Exchange diff level    : %lf %lf %lf %lf\n",
                   average[6+9*i], stddev[6+9*i], minimum[6+9*i],
                   maximum[6+9*i]);
            fprintf(fp, "     Apply BC               : %lf %lf %lf %lf\n",
                   average[7+9*i], stddev[7+9*i], minimum[7+9*i],
                   maximum[7+9*i]);
            fprintf(fp, "     Wait time              : %lf %lf %lf %lf\n",
                   average[8+9*i], stddev[8+9*i], minimum[8+9*i],
                   maximum[8+9*i]);
            fprintf(fp, "     Unpack faces           : %lf %lf %lf %lf\n\n",
                   average[9+9*i], stddev[9+9*i], minimum[9+9*i],
                   maximum[9+9*i]);

            if (!i) {
               fprintf(fp, "     Comm partners total ave: %lf %lf %lf %lf\n",
                       average[120], stddev[120], minimum[120], maximum[120]);
               fprintf(fp, "     Comm partners total min: %lf %lf %lf %lf\n",
                       average[125], stddev[125], minimum[125], maximum[125]);
               fprintf(fp, "     Comm partners total max: %lf %lf %lf %lf\n",
                       average[130], stddev[130], minimum[130], maximum[130]);
               fprintf(fp, "     Comm partners uniq ave : %lf %lf %lf %lf\n",
                       average[121], stddev[121], minimum[121], maximum[121]);
               fprintf(fp, "     Comm partners uniq min : %lf %lf %lf %lf\n",
                       average[126], stddev[126], minimum[126], maximum[126]);
               fprintf(fp, "     Comm partners uniq max : %lf %lf %lf %lf\n",
                       average[131], stddev[131], minimum[131], maximum[131]);
            } else {
               fprintf(fp, "     Comm partners average  : %lf %lf %lf %lf\n",
                       average[116+i], stddev[116+i], minimum[116+i],
                       maximum[116+i]);
               fprintf(fp, "     Comm partners minimum  : %lf %lf %lf %lf\n",
                       average[121+i], stddev[121+i], minimum[121+i],
                       maximum[121+i]);
               fprintf(fp, "     Comm partners maximum  : %lf %lf %lf %lf\n",
                       average[126+i], stddev[126+i], minimum[126+i],
                       maximum[126+i]);
            }
            fprintf(fp, "     Messages received      : %lf %lf %lf %lf\n",
               average[71+9*i], stddev[71+9*i], minimum[71+9*i],
                   maximum[71+9*i]);
            fprintf(fp, "     Bytes received         : %lf %lf %lf %lf\n",
               average[69+9*i], stddev[69+9*i], minimum[69+9*i],
                   maximum[69+9*i]);
            fprintf(fp, "     Faces received         : %lf %lf %lf %lf\n",
               average[73+9*i], stddev[73+9*i], minimum[73+9*i],
                   maximum[73+9*i]);
            fprintf(fp, "     Messages sent          : %lf %lf %lf %lf\n",
               average[72+9*i], stddev[72+9*i], minimum[72+9*i],
                   maximum[72+9*i]);
            fprintf(fp, "     Bytes sent             : %lf %lf %lf %lf\n",
               average[70+9*i], stddev[70+9*i], minimum[70+9*i],
                   maximum[70+9*i]);
            fprintf(fp, "     Faces sent             : %lf %lf %lf %lf\n",
               average[74+9*i], stddev[74+9*i], minimum[74+9*i],
                   maximum[74+9*i]);
            fprintf(fp, "     Faces exchanged same   : %lf %lf %lf %lf\n",
               average[76+9*i], stddev[76+9*i], minimum[76+9*i],
                   maximum[76+9*i]);
            fprintf(fp, "     Faces exchanged diff   : %lf %lf %lf %lf\n",
               average[77+9*i], stddev[77+9*i], minimum[77+9*i],
                   maximum[77+9*i]);
            fprintf(fp, "     Faces with BC applied  : %lf %lf %lf %lf\n",
               average[75+9*i], stddev[75+9*i], minimum[75+9*i],
                   maximum[75+9*i]);
         }

         fprintf(fp, "\n---------------------------------------------\n");
         fprintf(fp, "             Gridsum performance\n");
         fprintf(fp, "---------------------------------------------\n\n");
         fprintf(fp, "     Time: ave, stddev, min, max (sec): %lf %lf %lf %lf\n\n",
                average[39], stddev[39], minimum[39], maximum[39]);
         fprintf(fp, "        red : ave, stddev, min, max (sec): %lf %lf %lf %lf\n\n",
                average[40], stddev[40], minimum[40], maximum[40]);
         fprintf(fp, "        calc: ave, stddev, min, max (sec): %lf %lf %lf %lf\n\n",
                average[41], stddev[41], minimum[41], maximum[41]);
         fprintf(fp, "     total number:             %d\n", total_red);
         fprintf(fp, "     number per timestep:      %d\n\n", num_vars);

         fprintf(fp, "---------------------------------------------\n");
         fprintf(fp, "               Mesh Refinement\n");
         fprintf(fp, "---------------------------------------------\n\n");
         fprintf(fp, "     Time: ave, stddev, min, max (sec): %lf %lf %lf %lf\n\n",
                average[42], stddev[42], minimum[42], maximum[42]);
         fprintf(fp, "     Number of refinement steps: %d\n\n", nrs);
         fprintf(fp, "     Number of load balance steps: %d\n\n", nlbs);
         fprintf(fp, "     Number of redistributing steps: %d\n\n", nrrs);
         fprintf(fp, "     Total blocks           : %lld\n", total_blocks);
         fprintf(fp, "     Blocks/timestep ave, min, max : %lf %lld %lld\n",
                ((double) total_blocks)/((double) (num_tsteps*stages_per_ts)),
                (long long) nb_min, (long long) nb_max);
         fprintf(fp, "     Max blocks on a processor at any time: %d\n",
                global_max_b);
         fprintf(fp, "     total blocks split     : %lf\n",
                 average[105]*num_pes);
         fprintf(fp, "     total blocks reformed  : %lf\n\n",
                 average[106]*num_pes);
         fprintf(fp, "     total blocks moved     : %lf\n",
                 average[107]*num_pes);
         fprintf(fp, "     total moved load bal   : %lf\n",
                 average[108]*num_pes);
         fprintf(fp, "     total moved redistribut: %lf\n",
                 average[109]*num_pes);
         fprintf(fp, "     total moved coasening  : %lf\n\n",
                 average[110]*num_pes);
         fprintf(fp, "     parents ave, min, max  : %lf %d %d\n", average[132],
                 (int) minimum[132], (int) maximum[132]);
         fprintf(fp, "     max dots ave, min, max : %lf %d %d\n", average[133],
                 (int) minimum[133], (int) maximum[133]);
         fprintf(fp, "     ave dots used          : %lf\n\n",
                 average[134]/nlbs);
         fprintf(fp, "                              average    stddev  minimum  maximum\n");
         fprintf(fp, "     Per processor:\n");
         fprintf(fp, "     total blocks split     : %lf %lf %lf %lf\n",
                average[105], stddev[105], minimum[105], maximum[105]);
         fprintf(fp, "     total blocks reformed  : %lf %lf %lf %lf\n",
                average[106], stddev[106], minimum[106], maximum[106]);
         fprintf(fp, "     Total blocks moved     : %lf %lf %lf %lf\n",
                average[107], stddev[107], minimum[107], maximum[107]);
         fprintf(fp, "     Blocks moved load bal  : %lf %lf %lf %lf\n",
                average[108], stddev[108], minimum[108], maximum[108]);
         fprintf(fp, "     Blocks moved redistribu: %lf %lf %lf %lf\n",
                average[109], stddev[109], minimum[109], maximum[109]);
         fprintf(fp, "     Blocks moved coarsening: %lf %lf %lf %lf\n",
                average[110], stddev[110], minimum[110], maximum[110]);
         fprintf(fp, "     Time:\n");
         fprintf(fp, "        compare objects     : %lf %lf %lf %lf\n",
                average[43], stddev[43], minimum[43], maximum[43]);
         fprintf(fp, "        mark refine/coarsen : %lf %lf %lf %lf\n",
                average[44], stddev[44], minimum[44], maximum[44]);
         fprintf(fp, "        communicate block 1 : %lf %lf %lf %lf\n",
             average[47], stddev[47], minimum[47], maximum[47]);
         fprintf(fp, "        split blocks        : %lf %lf %lf %lf\n",
                average[46], stddev[46], minimum[46], maximum[46]);
         fprintf(fp, "        communicate block 2 : %lf %lf %lf %lf\n",
                average[48], stddev[48], minimum[48], maximum[48]);
         fprintf(fp, "        sync time           : %lf %lf %lf %lf\n",
                average[49], stddev[49], minimum[49], maximum[49]);
         fprintf(fp, "        misc time           : %lf %lf %lf %lf\n",
                average[45], stddev[45], minimum[45], maximum[45]);
         fprintf(fp, "        total coarsen blocks: %lf %lf %lf %lf\n",
                average[50], stddev[50], minimum[50], maximum[50]);
         fprintf(fp, "           coarsen blocks   : %lf %lf %lf %lf\n",
                average[51], stddev[51], minimum[51], maximum[51]);
         fprintf(fp, "           pack blocks      : %lf %lf %lf %lf\n",
                average[52], stddev[52], minimum[52], maximum[52]);
         fprintf(fp, "           move blocks      : %lf %lf %lf %lf\n",
                average[53], stddev[53], minimum[53], maximum[53]);
         fprintf(fp, "           unpack blocks    : %lf %lf %lf %lf\n",
                average[54], stddev[54], minimum[54], maximum[54]);
         fprintf(fp, "        total redistribute  : %lf %lf %lf %lf\n",
                average[64], stddev[64], minimum[64], maximum[64]);
         fprintf(fp, "           choose blocks    : %lf %lf %lf %lf\n",
                average[65], stddev[65], minimum[65], maximum[65]);
         fprintf(fp, "           pack blocks      : %lf %lf %lf %lf\n",
                average[66], stddev[66], minimum[66], maximum[66]);
         fprintf(fp, "           move blocks      : %lf %lf %lf %lf\n",
                average[67], stddev[67], minimum[67], maximum[67]);
         fprintf(fp, "           unpack blocks    : %lf %lf %lf %lf\n",
                average[68], stddev[68], minimum[68], maximum[68]);
         fprintf(fp, "        total load balance  : %lf %lf %lf %lf\n",
                average[55], stddev[55], minimum[55], maximum[55]);
         fprintf(fp, "           sort             : %lf %lf %lf %lf\n",
                average[56], stddev[56], minimum[56], maximum[56]);
         fprintf(fp, "           move dots back   : %lf %lf %lf %lf\n",
                average[61], stddev[61], minimum[61], maximum[61]);
         fprintf(fp, "           move blocks total: %lf %lf %lf %lf\n",
                average[62], stddev[62], minimum[62], maximum[62]);
         fprintf(fp, "              pack blocks   : %lf %lf %lf %lf\n",
                average[57], stddev[57], minimum[57], maximum[57]);
         fprintf(fp, "              move blocks   : %lf %lf %lf %lf\n",
                average[58], stddev[58], minimum[58], maximum[58]);
         fprintf(fp, "              unpack blocks : %lf %lf %lf %lf\n",
                average[59], stddev[59], minimum[59], maximum[59]);
         fprintf(fp, "              misc          : %lf %lf %lf %lf\n\n",
                average[60], stddev[60], minimum[60], maximum[60]);

         fprintf(fp, "---------------------------------------------\n");
         fprintf(fp, "                   Plot\n");
         fprintf(fp, "---------------------------------------------\n\n");
         fprintf(fp, "     Time: ave, stddev, min, max (sec): %lf %lf %lf %lf\n\n",
                average[63], stddev[63], minimum[63], maximum[63]);
         fprintf(fp, "     Number of plot steps: %d\n", nps);
         fprintf(fp, "\n ================== End report ===================\n");

         fclose(fp);
      }

      if (report_perf & 4) {
         printf("\n ================ Start report ===================\n\n");
         printf("          Mantevo miniAMR\n");
         printf("          version %s\n\n", version);

         printf("Run on %d ranks arranged in a %d x %d x %d grid\n", num_pes,
                npx, npy, npz);
         printf("initial blocks per rank %d x %d x %d\n", init_block_x,
                init_block_y, init_block_z);
         printf("block size %d x %d x %d\n", x_block_size, y_block_size,
                z_block_size);
         if (reorder)
            printf("Initial ranks arranged by RCB across machine\n\n");
         else
            printf("Initial ranks arranged as a grid across machine\n\n");
         printf("change_dir %d group_blocks %d limit_move %d\n", change_dir,
                group_blocks, limit_move);
         if (permute)
            printf("Order of exchanges permuted\n");
         printf("Maximum number of blocks per rank is %d\n", max_num_blocks);
         if (code)
            printf("Code set to code %d\n", code);
         printf("Number of levels of refinement is %d\n", num_refine);
         printf("Blocks can change by %d levels per refinement step\n",
            block_change);
         if (refine_ghost)
            printf("Ghost cells will be used determine is block is refined\n");
         if (uniform_refine)
            printf("\nBlocks will be uniformly refined\n");
         else {
            printf("\nBlocks will be refined by %d objects\n\n", num_objects);
            for (i = 0; i < num_objects; i++) {
               op = &objects[i];
               if (op->type == 0)
                  printf("Object %d is the surface of a rectangle\n", i);
               else if (op->type == 1)
                  printf("Object %d is the volume of a rectangle\n", i);
               else if (op->type == 2)
                  printf("Object %d is the surface of a spheroid\n", i);
               else if (op->type == 3)
                  printf("Object %d is the volume of a spheroid\n", i);
               else if (op->type == 4)
                  printf("Object %d is the surface of x+ hemispheroid\n", i);
               else if (op->type == 5)
                  printf("Object %d is the volume of x+ hemispheroid\n", i);
               else if (op->type == 6)
                  printf("Object %d is the surface of x- hemispheroid\n", i);
               else if (op->type == 7)
                  printf("Object %d is the volume of x- hemispheroid\n", i);
               else if (op->type == 8)
                  printf("Object %d is the surface of y+ hemispheroid\n", i);
               else if (op->type == 9)
                  printf("Object %d is the volume of y+ hemispheroid\n", i);
               else if (op->type == 10)
                  printf("Object %d is the surface of y- hemispheroid\n", i);
               else if (op->type == 11)
                  printf("Object %d is the volume of y- hemispheroid\n", i);
               else if (op->type == 12)
                  printf("Object %d is the surface of z+ hemispheroid\n", i);
               else if (op->type == 13)
                  printf("Object %d is the volume of z+ hemispheroid\n", i);
               else if (op->type == 14)
                  printf("Object %d is the surface of z- hemispheroid\n", i);
               else if (op->type == 15)
                  printf("Object %d is the volume of z- hemispheroid\n", i);
               else if (op->type == 20)
                  printf("Object %d is the surface of x axis cylinder\n", i);
               else if (op->type == 21)
                  printf("Object %d is the volune of x axis cylinder\n", i);
               else if (op->type == 22)
                  printf("Object %d is the surface of y axis cylinder\n", i);
               else if (op->type == 23)
                  printf("Object %d is the volune of y axis cylinder\n", i);
               else if (op->type == 24)
                  printf("Object %d is the surface of z axis cylinder\n", i);
               else if (op->type == 25)
                  printf("Object %d is the volune of z axis cylinder\n", i);
               if (op->bounce == 0)
                  printf("Object may leave mesh\n");
               else
                  printf("Object center will bounce off of walls\n");
               printf("Center starting at %lf %lf %lf\n",
                      op->orig_cen[0], op->orig_cen[1], op->orig_cen[2]);
               printf("Center end at %lf %lf %lf\n",
                      op->cen[0], op->cen[1], op->cen[2]);
               if (use_time) {
                  printf("Velocity of %lf %lf %lf\n",
                         op->orig_move[0], op->orig_move[1], op->orig_move[2]);
                  delta = end_time/((double) num_tsteps);
                  printf("   Rate relative to smallest cell size %lf %lf %lf\n",
                         delta*op->orig_move[0]*((double)
                                                 (mesh_size[0]*x_block_size)),
                         delta*op->orig_move[1]*((double)
                                                 (mesh_size[1]*y_block_size)),
                         delta*op->orig_move[2]*((double)
                                                 (mesh_size[2]*z_block_size)));
               } else {
                  printf("Moving at %lf %lf %lf per timestep\n",
                         op->orig_move[0], op->orig_move[1], op->orig_move[2]);
                  printf("   Rate relative to smallest cell size %lf %lf %lf\n",
                      op->orig_move[0]*((double) (mesh_size[0]*x_block_size)),
                      op->orig_move[1]*((double) (mesh_size[1]*y_block_size)),
                      op->orig_move[2]*((double) (mesh_size[2]*z_block_size)));
               }
               printf("Initial size %lf %lf %lf\n",
                      op->orig_size[0], op->orig_size[1], op->orig_size[2]);
               printf("Final size %lf %lf %lf\n",
                      op->size[0], op->size[1], op->size[2]);
               printf("Size increasing %lf %lf %lf per timestep\n",
                      op->inc[0], op->inc[1], op->inc[2]);
               printf("   Rate relative to smallest cell size %lf %lf %lf\n\n",
                      op->inc[0]*((double) (mesh_size[0]*x_block_size)),
                      op->inc[1]*((double) (mesh_size[1]*y_block_size)),
                      op->inc[2]*((double) (mesh_size[2]*z_block_size)));
            }
         }
         if (use_time)
            printf("\nTime %lf in %d timesteps\n", end_time, num_tsteps);
         else
            printf("\nNumber of timesteps is %d\n", num_tsteps);
         printf("Communicaion/computation stages per timestep is %d\n",
                stages_per_ts);
         if (nonblocking)
            printf("Communication will be performed with nonblocking sends\n");
         else
            printf("Communication will be performed with blocking sends\n");
         printf("Will perform checksums every %d stages\n", checksum_freq);
         printf("Will refine every %d timesteps\n", refine_freq);
         if (!lb_method)
            printf("Load balance by RCB (Recursive Coordinate Bisection)\n");
         else if (lb_method == 1)
            printf("Load balance by Morton Space Filling Curve\n");
         else if (lb_method == 2)
            printf("Load balance by Peano Hilbert Space Filling Curve\n");
         else
            printf("Load balance by truncated Hilbert style Space Filling Curve\n");
         if (lb_opt == 0)
            printf("Load balance will not be performed\n");
         else
            printf("Load balance when inbalanced by %d%%\n", inbalance);
         if (lb_opt == 2)
            printf("Load balance at each phase of refinement step\n");
         if (plot_freq)
            printf("Will plot results every %d timesteps\n", plot_freq);
         else
            printf("Will not plot results\n");
         if (stencil)
            printf("Calculate on %d variables with %d point stencil\n",
                   num_vars, stencil);
         else
            printf("Calculate on %d variables with variable stencils\n",
                   num_vars);
         printf("Communicate %d variables at a time\n", comm_vars);
         printf("Error tolorance for variable sums is 10^(-%d)\n", error_tol);
         if (send_faces)
            printf("Will send data from each face seperately\n");

         printf("\nTotal time for test: ave, std, min, max (sec): %lf %lf %lf %lf\n\n",
                average[0], stddev[0], minimum[0], maximum[0]);

         printf("\nNumber of malloc calls: ave, std, min, max (sec): %lf %lf %lf %lf\n",
                average[111], stddev[111], minimum[111], maximum[111]);
         printf("\nAmount malloced: ave, std, min, max: %lf %lf %lf %lf\n",
                average[112], stddev[112], minimum[112], maximum[112]);
         printf("\nMalloc calls in init: ave, std, min, max (sec): %lf %lf %lf %lf\n",
                average[113], stddev[113], minimum[113], maximum[113]);
         printf("\nAmount malloced in init: ave, std, min, max: %lf %lf %lf %lf\n",
                average[114], stddev[114], minimum[114], maximum[114]);
         printf("\nMalloc calls in timestepping: ave, std, min, max (sec): %lf %lf %lf %lf\n",
                average[115], stddev[115], minimum[115], maximum[115]);
         printf("\nAmount malloced in timestepping: ave, std, min, max: %lf %lf %lf %lf\n\n",
                average[116], stddev[116], minimum[116], maximum[116]);

         printf("Main (parse, allocate) Time: ave, stddev, min, max (sec): %lf %lf %lf %lf\n\n",
                average[138], stddev[138], minimum[138], maximum[138]);
         printf("Initailize Time: ave, stddev, min, max (sec): %lf %lf %lf %lf\n\n",
                average[137], stddev[137], minimum[137], maximum[137]);
         printf("---------------------------------------------\n");
         printf("          Computational Performance\n");
         printf("---------------------------------------------\n\n");
         printf("     Time: ave, stddev, min, max (sec): %lf %lf %lf %lf\n\n",
                average[38], stddev[38], minimum[38], maximum[38]);
         printf("     total GFLOPS:             %lf\n", total_gflops);
         printf("     Average GFLOPS per rank:  %lf\n\n", gflops_rank);
         printf("     Total floating point ops: %lf\n\n", total_fp_ops);
         printf("        Adds:                  %lf\n", average[139]);
         printf("        Muls:                  %lf\n", average[140]);
         printf("        Divides:               %lf\n\n", average[141]);
         printf("     Sum of min/max compute times per ts: %lf %lf\n\n",
                tmin[1], tmax[1]);

         printf("---------------------------------------------\n");
         printf("           Interblock communication\n");
         printf("---------------------------------------------\n\n");
         printf("     Time: ave, stddev, min, max (sec): %lf %lf %lf %lf\n\n",
                average[37], stddev[37], minimum[37], maximum[37]);
         for (i = 0; i < 4; i++) {
            if (i == 0)
               printf("\nTotal communication:\n\n");
            else if (i == 1)
               printf("\nX direction communication statistics:\n\n");
            else if (i == 2)
               printf("\nY direction communication statistics:\n\n");
            else
               printf("\nZ direction communication statistics:\n\n");
            printf("                              average    stddev  minimum  maximum\n");
            printf("     Total                  : %lf %lf %lf %lf\n",
                   average[1+9*i], stddev[1+9*i], minimum[1+9*i],
                   maximum[1+9*i]);
            printf("     Post IRecv             : %lf %lf %lf %lf\n",
                   average[2+9*i], stddev[2+9*i], minimum[2+9*i],
                   maximum[2+9*i]);
            printf("     Pack faces             : %lf %lf %lf %lf\n",
                   average[3+9*i], stddev[3+9*i], minimum[3+9*i],
                   maximum[3+9*i]);
            printf("     Send messages          : %lf %lf %lf %lf\n",
                   average[4+9*i], stddev[4+9*i], minimum[4+9*i],
                   maximum[4+9*i]);
            printf("     Exchange same level    : %lf %lf %lf %lf\n",
                   average[5+9*i], stddev[5+9*i], minimum[5+9*i],
                   maximum[5+9*i]);
            printf("     Exchange diff level    : %lf %lf %lf %lf\n",
                   average[6+9*i], stddev[6+9*i], minimum[6+9*i],
                   maximum[6+9*i]);
            printf("     Apply BC               : %lf %lf %lf %lf\n",
                   average[7+9*i], stddev[7+9*i], minimum[7+9*i],
                   maximum[7+9*i]);
            printf("     Wait time              : %lf %lf %lf %lf\n",
                   average[8+9*i], stddev[8+9*i], minimum[8+9*i],
                   maximum[8+9*i]);
            printf("     Unpack faces           : %lf %lf %lf %lf\n\n",
                   average[9+9*i], stddev[9+9*i], minimum[9+9*i],
                   maximum[9+9*i]);

            if (!i) {
               printf("     Comm partners total ave: %lf %lf %lf %lf\n",
                      average[120], stddev[120], minimum[120], maximum[120]);
               printf("     Comm partners total min: %lf %lf %lf %lf\n",
                      average[125], stddev[125], minimum[125], maximum[125]);
               printf("     Comm partners total max: %lf %lf %lf %lf\n",
                      average[130], stddev[130], minimum[130], maximum[130]);
               printf("     Comm partners uniq ave : %lf %lf %lf %lf\n",
                      average[121], stddev[121], minimum[121], maximum[121]);
               printf("     Comm partners uniq min : %lf %lf %lf %lf\n",
                      average[126], stddev[126], minimum[126], maximum[126]);
               printf("     Comm partners uniq max : %lf %lf %lf %lf\n",
                      average[131], stddev[131], minimum[131], maximum[131]);
            } else {
               printf("     Comm partners average  : %lf %lf %lf %lf\n",
                      average[116+i], stddev[116+i], minimum[116+i],
                      maximum[116+i]);
               printf("     Comm partners minimum  : %lf %lf %lf %lf\n",
                      average[121+i], stddev[121+i], minimum[121+i],
                      maximum[121+i]);
               printf("     Comm partners maximum  : %lf %lf %lf %lf\n",
                      average[126+i], stddev[126+i], minimum[126+i],
                      maximum[126+i]);
            }
            printf("     Messages received      : %lf %lf %lf %lf\n",
                   average[71+9*i], stddev[71+9*i], minimum[71+9*i],
                   maximum[71+9*i]);
            printf("     Bytes received         : %lf %lf %lf %lf\n",
                   average[69+9*i], stddev[69+9*i], minimum[69+9*i],
                   maximum[69+9*i]);
            printf("     Faces received         : %lf %lf %lf %lf\n",
                   average[73+9*i], stddev[73+9*i], minimum[73+9*i],
                   maximum[73+9*i]);
            printf("     Messages sent          : %lf %lf %lf %lf\n",
                   average[72+9*i], stddev[72+9*i], minimum[72+9*i],
                   maximum[72+9*i]);
            printf("     Bytes sent             : %lf %lf %lf %lf\n",
                   average[70+9*i], stddev[70+9*i], minimum[70+9*i],
                   maximum[70+9*i]);
            printf("     Faces sent             : %lf %lf %lf %lf\n",
                   average[74+9*i], stddev[74+9*i], minimum[74+9*i],
                   maximum[74+9*i]);
            printf("     Faces exchanged same   : %lf %lf %lf %lf\n",
                   average[76+9*i], stddev[76+9*i], minimum[76+9*i],
                   maximum[76+9*i]);
            printf("     Faces exchanged diff   : %lf %lf %lf %lf\n",
                   average[77+9*i], stddev[77+9*i], minimum[77+9*i],
                   maximum[77+9*i]);
            printf("     Faces with BC applied  : %lf %lf %lf %lf\n",
                   average[75+9*i], stddev[75+9*i], minimum[75+9*i],
                   maximum[75+9*i]);
         }
         printf("\n     Sum of min/max communicate times per ts: %lf %lf\n",
                tmin[0], tmax[0]);
         printf("\n     Sum of min/max calc and comm times per ts: %lf %lf\n",
                tmin[4], tmax[4]);

         printf("\n---------------------------------------------\n");
         printf("             Gridsum performance\n");
         printf("---------------------------------------------\n\n");
         printf("     Time: ave, stddev, min, max (sec): %lf %lf %lf %lf\n\n",
                average[39], stddev[39], minimum[39], maximum[39]);
         printf("        red : ave, stddev, min, max (sec): %lf %lf %lf %lf\n\n",
                average[40], stddev[40], minimum[40], maximum[40]);
         printf("        calc: ave, stddev, min, max (sec): %lf %lf %lf %lf\n\n",
                average[41], stddev[41], minimum[41], maximum[41]);
         printf("     total number:             %d\n", total_red);
         printf("     number per timestep:      %d\n\n", num_vars);

         printf("     Sum of min/max gridsum times per ts: %lf %lf\n\n",
                tmin[2], tmax[2]);

         printf("---------------------------------------------\n");
         printf("               Mesh Refinement\n");
         printf("---------------------------------------------\n\n");
         printf("     Time: ave, stddev, min, max (sec): %lf %lf %lf %lf\n\n",
                average[42], stddev[42], minimum[42], maximum[42]);
         printf("     Number of refinement steps: %d\n\n", nrs);
         printf("     Number of load balance steps: %d\n\n", nlbs);
         printf("     Number of redistributing steps: %d\n\n", nrrs);
         printf("     Total blocks           : %lld\n", total_blocks);
         printf("     Blocks/timestep ave, min, max : %lf %lld %lld\n",
                ((double) total_blocks)/((double) (num_tsteps*stages_per_ts)),
                (long long) nb_min, (long long) nb_max);
         printf("     Max blocks on a processor at any time: %d\n",
                global_max_b);
         printf("     total blocks split     : %lf\n", average[105]*num_pes);
         printf("     total blocks reformed  : %lf\n\n", average[106]*num_pes);
         printf("     total blocks moved     : %lf\n", average[107]*num_pes);
         printf("     total moved load bal   : %lf\n", average[108]*num_pes);
         printf("     total moved redistribut: %lf\n", average[109]*num_pes);
         printf("     total moved coasening  : %lf\n\n", average[110]*num_pes);
         printf("     parents ave, min, max  : %lf %d %d\n", average[132],
                (int) minimum[132], (int) maximum[132]);
         printf("     max dots ave, min, max : %lf %d %d\n", average[133],
                (int) minimum[133], (int) maximum[133]);
         printf("     ave dots used          : %lf\n\n", average[134]/nlbs);
         printf("     number of times more sets of blocks than ranks: %d\n",
                num_over);
         printf("     total sets of blocks over num of ranks: %d\n", tot_over);
         printf("     max number of sets of blocks on a rank: %d\n\n",
                max_groups);
         printf("                              average    stddev  minimum  maximum\n");
         printf("     Per processor:\n");
         printf("     total blocks split     : %lf %lf %lf %lf\n",
                average[105], stddev[105], minimum[105], maximum[105]);
         printf("     total blocks reformed  : %lf %lf %lf %lf\n",
                average[106], stddev[106], minimum[106], maximum[106]);
         printf("     Total blocks moved     : %lf %lf %lf %lf\n",
                average[107], stddev[107], minimum[107], maximum[107]);
         printf("     Blocks moved load bal  : %lf %lf %lf %lf\n",
                average[108], stddev[108], minimum[108], maximum[108]);
         printf("     Blocks moved redistribu: %lf %lf %lf %lf\n",
                average[109], stddev[109], minimum[109], maximum[109]);
         printf("     Blocks moved coarsening: %lf %lf %lf %lf\n",
                average[110], stddev[110], minimum[110], maximum[110]);
         printf("     Time:\n");
         printf("        initial refine      : %lf %lf %lf %lf\n\n",
                average[136], stddev[136], minimum[136], maximum[136]);
         printf("        compare objects     : %lf %lf %lf %lf\n",
                average[43], stddev[43], minimum[43], maximum[43]);
         printf("        mark refine/coarsen : %lf %lf %lf %lf\n",
                average[44], stddev[44], minimum[44], maximum[44]);
         printf("        communicate block 1 : %lf %lf %lf %lf\n",
             average[47], stddev[47], minimum[47], maximum[47]);
         printf("        split blocks        : %lf %lf %lf %lf\n",
                average[46], stddev[46], minimum[46], maximum[46]);
         printf("        communicate block 2 : %lf %lf %lf %lf\n",
                average[48], stddev[48], minimum[48], maximum[48]);
         printf("        sync time           : %lf %lf %lf %lf\n",
                average[49], stddev[49], minimum[49], maximum[49]);
         printf("        time for groups     : %lf %lf %lf %lf\n",
                average[135], stddev[135], minimum[135], maximum[135]);
         printf("        misc time           : %lf %lf %lf %lf\n",
                average[45], stddev[45], minimum[45], maximum[45]);
         printf("        total coarsen blocks: %lf %lf %lf %lf\n",
                average[50], stddev[50], minimum[50], maximum[50]);
         printf("           coarsen blocks   : %lf %lf %lf %lf\n",
                average[51], stddev[51], minimum[51], maximum[51]);
         printf("           pack blocks      : %lf %lf %lf %lf\n",
                average[52], stddev[52], minimum[52], maximum[52]);
         printf("           move blocks      : %lf %lf %lf %lf\n",
                average[53], stddev[53], minimum[53], maximum[53]);
         printf("           unpack blocks    : %lf %lf %lf %lf\n",
                average[54], stddev[54], minimum[54], maximum[54]);
         printf("        total redistribute  : %lf %lf %lf %lf\n",
                average[64], stddev[64], minimum[64], maximum[64]);
         printf("           choose blocks    : %lf %lf %lf %lf\n",
                average[65], stddev[65], minimum[65], maximum[65]);
         printf("           pack blocks      : %lf %lf %lf %lf\n",
                average[66], stddev[66], minimum[66], maximum[66]);
         printf("           move blocks      : %lf %lf %lf %lf\n",
                average[67], stddev[67], minimum[67], maximum[67]);
         printf("           unpack blocks    : %lf %lf %lf %lf\n",
                average[68], stddev[68], minimum[68], maximum[68]);
         printf("        total load balance  : %lf %lf %lf %lf\n",
                average[55], stddev[55], minimum[55], maximum[55]);
         printf("           sort             : %lf %lf %lf %lf\n",
                average[56], stddev[56], minimum[56], maximum[56]);
         printf("           move dots back   : %lf %lf %lf %lf\n",
                average[61], stddev[61], minimum[61], maximum[61]);
         printf("           move blocks total: %lf %lf %lf %lf\n",
                average[62], stddev[62], minimum[62], maximum[62]);
         printf("              pack blocks   : %lf %lf %lf %lf\n",
                average[57], stddev[57], minimum[57], maximum[57]);
         printf("              move blocks   : %lf %lf %lf %lf\n",
                average[58], stddev[58], minimum[58], maximum[58]);
         printf("              unpack blocks : %lf %lf %lf %lf\n",
                average[59], stddev[59], minimum[59], maximum[59]);
         printf("              misc          : %lf %lf %lf %lf\n\n",
                average[60], stddev[60], minimum[60], maximum[60]);

         printf("     Sum of min/max refinement times per ts: %lf %lf\n\n",
                tmin[3], tmax[3]);

         printf("---------------------------------------------\n");
         printf("                   Plot\n");
         printf("---------------------------------------------\n\n");
         printf("     Time: ave, stddev, min, max (sec): %lf %lf %lf %lf\n\n",
                average[63], stddev[63], minimum[63], maximum[63]);
         printf("     Number of plot steps: %d\n", nps);
         printf("\n ================== End report ===================\n");
printf("Summary: ranks %d ts %d time %lf calc %lf comm %lf red %lf refine %lf blocks/ts %lf max_blocks %d\n", num_pes, num_tsteps, average[0], average[38], average[37], average[39], average[42], ((double) total_blocks)/((double) (num_tsteps*stages_per_ts)), global_max_b);
fflush(NULL);
      }
   }
}

void calculate_results(void)
{
   double results[142], stddev_sum[139];
   int i;

   MPI_Allreduce(&local_max_b, &global_max_b, 1, MPI_INT, MPI_MAX,
                 MPI_COMM_WORLD);
   results[0] = timer_all;
   for (i = 0; i < 9; i++)
      results[i+1] = 0.0;
   for (i = 0; i < 3; i++) {
      results[1] += results[10+9*i] = timer_comm_dir[i];
      results[2] += results[11+9*i] = timer_comm_recv[i];
      results[3] += results[12+9*i] = timer_comm_pack[i];
      results[4] += results[13+9*i] = timer_comm_send[i];
      results[5] += results[14+9*i] = timer_comm_same[i];
      results[6] += results[15+9*i] = timer_comm_diff[i];
      results[7] += results[16+9*i] = timer_comm_bc[i];
      results[8] += results[17+9*i] = timer_comm_wait[i];
      results[9] += results[18+9*i] = timer_comm_unpack[i];
   }
   results[37] = timer_comm_all;
   results[38] = timer_calc_all;
   results[39] = timer_cs_all;
   results[40] = timer_cs_red;
   results[41] = timer_cs_calc;
   results[42] = timer_refine_all;
   results[43] = timer_refine_co;
   results[44] = timer_refine_mr;
   results[45] = timer_refine_cc;
   results[46] = timer_refine_sb;
   results[47] = timer_refine_c1;
   results[48] = timer_refine_c2;
   results[49] = timer_refine_sy;
   results[50] = timer_cb_all;
   results[51] = timer_cb_cb;
   results[52] = timer_cb_pa;
   results[53] = timer_cb_mv;
   results[54] = timer_cb_un;
   results[55] = timer_lb_all;
   results[56] = timer_lb_sort;
   results[57] = timer_lb_pa;
   results[58] = timer_lb_mv;
   results[59] = timer_lb_un;
   results[60] = timer_lb_misc;
   results[61] = timer_lb_mb;
   results[62] = timer_lb_ma;
   results[63] = timer_plot;
   results[64] = timer_rs_all;
   results[65] = timer_rs_ca;
   results[66] = timer_rs_pa;
   results[67] = timer_rs_mv;
   results[68] = timer_rs_un;
   for (i = 0; i < 9; i++)
      results[69+i] = 0.0;
   for (i = 0; i < 3; i++) {
      results[69] += results[78+9*i] = size_mesg_recv[i];
      results[70] += results[79+9*i] = size_mesg_send[i];
      results[71] += results[80+9*i] = (double) counter_halo_recv[i];
      results[72] += results[81+9*i] = (double) counter_halo_send[i];
      results[73] += results[82+9*i] = (double) counter_face_recv[i];
      results[74] += results[83+9*i] = (double) counter_face_send[i];
      results[75] += results[84+9*i] = (double) counter_bc[i];
      results[76] += results[85+9*i] = (double) counter_same[i];
      results[77] += results[86+9*i] = (double) counter_diff[i];
   }
   results[105] = (double) num_refined;
   results[106] = (double) num_reformed;
   num_moved_all = num_moved_lb + num_moved_coarsen + num_moved_rs;
   results[107] = (double) num_moved_all;
   results[108] = (double) num_moved_lb;
   results[109] = (double) num_moved_rs;
   results[110] = (double) num_moved_coarsen;
   results[111] = (double) counter_malloc;
   results[112] = size_malloc;
   results[113] = (double) counter_malloc_init;
   results[114] = size_malloc_init;
   results[115] = (double) (counter_malloc - counter_malloc_init);
   results[116] = size_malloc - size_malloc_init;
   results[117] = (double) num_comm_x/(double) nrs;
   results[118] = (double) num_comm_y/(double) nrs;
   results[119] = (double) num_comm_z/(double) nrs;
   results[120] = (double) num_comm_tot/(double) nrs;
   results[121] = (double) num_comm_uniq/(double) nrs;
   results[122] = (double) num_comm_x_min;
   results[123] = (double) num_comm_y_min;
   results[124] = (double) num_comm_z_min;
   results[125] = (double) num_comm_t_min;
   results[126] = (double) num_comm_u_min;
   results[127] = (double) num_comm_x_max;
   results[128] = (double) num_comm_y_max;
   results[129] = (double) num_comm_z_max;
   results[130] = (double) num_comm_t_max;
   results[131] = (double) num_comm_u_max;
   results[132] = (double) num_parents;
   results[133] = (double) max_dots_used;
   results[134] = (double) total_dots_used;
   results[135] = (double) timer_group;
   results[136] = timer_refine_init;
   results[137] = timer_init;
   results[138] = timer_main;
   results[139] = total_fp_adds;
   results[140] = total_fp_muls;
   results[141] = total_fp_divs;

   MPI_Allreduce(results, average, 142, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(results, minimum, 139, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
   MPI_Allreduce(results, maximum, 139, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

   for (i = 0; i < 139; i++) {
      average[i] /= (double) num_pes;
      stddev[i] = (results[i] - average[i])*(results[i] - average[i]);
   }
   MPI_Allreduce(stddev, stddev_sum, 139, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   for (i = 0; i < 139; i++)
      stddev[i] = sqrt(stddev_sum[i]/((double) num_pes));
}

void init_profile(void)
{
   int i;

   timer_all = 0.0;

   timer_comm_all = 0.0;
   for (i = 0; i < 3; i++) {
      timer_comm_dir[i] = 0.0;
      timer_comm_recv[i] = 0.0;
      timer_comm_pack[i] = 0.0;
      timer_comm_send[i] = 0.0;
      timer_comm_same[i] = 0.0;
      timer_comm_diff[i] = 0.0;
      timer_comm_bc[i] = 0.0;
      timer_comm_wait[i] = 0.0;
      timer_comm_unpack[i] = 0.0;
   }

   timer_calc_all = 0.0;

   timer_cs_all = 0.0;
   timer_cs_red = 0.0;
   timer_cs_calc = 0.0;

   timer_refine_all = 0.0;
   timer_refine_co = 0.0;
   timer_refine_mr = 0.0;
   timer_refine_cc = 0.0;
   timer_refine_sb = 0.0;
   timer_refine_c1 = 0.0;
   timer_refine_c2 = 0.0;
   timer_refine_sy = 0.0;
   timer_refine_init = 0.0;
   timer_cb_all = 0.0;
   timer_cb_cb = 0.0;
   timer_cb_pa = 0.0;
   timer_cb_mv = 0.0;
   timer_cb_un = 0.0;
   timer_lb_all = 0.0;
   timer_lb_sort = 0.0;
   timer_lb_pa = 0.0;
   timer_lb_mv = 0.0;
   timer_lb_un = 0.0;
   timer_lb_misc = 0.0;
   timer_lb_mb = 0.0;
   timer_lb_ma = 0.0;
   timer_rs_all = 0.0;
   timer_rs_ca = 0.0;
   timer_rs_pa = 0.0;
   timer_rs_mv = 0.0;
   timer_rs_un = 0.0;

   timer_plot = 0.0;
   timer_init = 0.0;

   total_blocks = 0;
   max_dots_used = 0;
   total_dots_used = 0;
   nrrs = 0;
   nrs = 0;
   nps = 0;
   nlbs = 0;
   num_refined = 0;
   num_reformed = 0;
   num_moved_all = 0;
   num_moved_lb = 0;
   num_moved_rs = 0;
   num_moved_coarsen = 0;
   num_comm_x = 0;
   num_comm_y = 0;
   num_comm_z = 0;
   num_comm_tot = 0;
   num_comm_uniq = 0;
   num_comm_x_max = 0;
   num_comm_y_max = 0;
   num_comm_z_max = 0;
   num_comm_t_max = 0;
   num_comm_u_max = 0;
   num_comm_x_min = num_pes;
   num_comm_y_min = num_pes;
   num_comm_z_min = num_pes;
   num_comm_t_min = num_pes;
   num_comm_u_min = num_pes;
   for (i = 0; i < 3; i++) {
      counter_halo_recv[i] = 0;
      counter_halo_send[i] = 0;
      size_mesg_recv[i] = 0.0;
      size_mesg_send[i] = 0.0;
      counter_face_recv[i] = 0;
      counter_face_send[i] = 0;
      counter_bc[i] = 0;
      counter_same[i] = 0;
      counter_diff[i] = 0;
   }
   total_red = 0;
   num_over = 0;
   tot_over = 0;
   max_groups = 1;
   for (i = 1; i < 5; i++)
      tmin[i] = tmax[i] = 0.0;
}
