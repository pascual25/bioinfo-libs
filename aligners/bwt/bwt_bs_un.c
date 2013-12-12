#include "bwt.h"

//-----------------------------------------------------------------------------

#define jump_error      5

size_t bwt_map_inexact_read_bs_un(fastq_read_t *read, 
				  bwt_optarg_t *bwt_optarg, 
				  bwt_index_t *index, 
				  array_list_t *mapping_list_tmp,
				  int type) {

     alignment_t *alignment;
     char *seq = read->sequence;
     size_t len = read->length;

     const int MAX_BWT_ALIGNMENTS = 10;
     int filter_exceeded = 0;
     int pos;
     size_t num_mappings = 0;
     size_t idx, key, direction;
     size_t cigar_len, num_cigar_ops;
     size_t best_pos, array_size;
     size_t *allocate_pos_alignments;
     size_t k_start, l_start;
     size_t start_mapping;
     size_t tot_alignments = 0;
     char plusminus[2] = "-+";
     char error;
     results_list r_list;
     result *r;
     char *cigar_dup;
     char error_debug = 0;
     char cigar[1024];
     array_list_t *mapping_list_filter;
     alignment_t *best_alignment, *aux_alignment;

     char *seq_dup, *seq_strand;
     size_t start = 0;
     size_t end = len - 1;
     size_t len_calc = len;



     if (len < 5) {
	  char aux[len + 2];
	  sprintf(aux, "%luX", len);

	  char *quality_clipping = (char *)malloc(sizeof(char)*50);
	  sprintf(quality_clipping, "%i", NONE_HARD_CLIPPING);

	  alignment = alignment_new();
	  alignment_init_single_end(NULL,
				    strdup(seq),
				    quality_clipping,
				    0,
				    -1,
				    -1,
				    strdup(aux), 1, 0, 0, 0, 0, NULL, alignment);
	  array_list_insert((void*) alignment, mapping_list_tmp);

	  return 1;
     }
     
     char *code_seq = (char *) malloc(len * sizeof(char));

     bwt_encode_Bases(code_seq, seq, len, &index->table);

     // calculate vectors k and l
     size_t *k1 = (size_t *) malloc(len * sizeof(size_t));
     size_t *l1 = (size_t *) malloc(len * sizeof(size_t));
     size_t *ki1 = (size_t *) malloc(len * sizeof(size_t));
     size_t *li1 = (size_t *) malloc(len * sizeof(size_t));

     size_t last_k1, last_l1;
     size_t last_ki1, last_li1;

     int back1_nt, forw1_nt;

     seq_strand = strdup(seq);
     error = MISMATCH;
     
     char *quality_clipping = (char *) malloc(sizeof(char) * 50);
     seq_dup = (char *) malloc(sizeof(char) * (len + 1));
     
     new_results_list(&r_list, bwt_optarg->filter_read_mappings);
     
     array_list_t *tmp_mapping_list  = array_list_new(bwt_optarg->filter_read_mappings, 1.25f,
						      COLLECTION_MODE_ASYNCHRONIZED);
     
     array_list_t *mapping_list = array_list_new(bwt_optarg->filter_read_mappings, 1.25f,
						 COLLECTION_MODE_ASYNCHRONIZED);
     array_list_t *anchor_list = array_list_new(bwt_optarg->filter_read_mappings, 1.25f,
						COLLECTION_MODE_ASYNCHRONIZED);
     cal_t *cal;
     bwt_anchor_t * anchor;
     
     int min_size = 20;
     size_t first_loop = 0;
     ///////////////////////////////
     while ((int)end - (int)start >= min_size) {
     ///////////////////////////////
       //printf("Start - end (%lu : %lu)\n", start, end);
       //printf("gap %i (min gap %i)\n", (int)end - (int)start, min_size);
       LOG_DEBUG_F("Iteracion %lu, (%lu : %lu) from %s\n", first_loop, start, end, seq);
       BWExactSearchVectorBackward(code_seq, start, end, 0, index->h_O.siz - 2,
				   k1, l1, &index->h_C, &index->h_C1, &index->h_O,
				   &last_k1, &last_l1, &back1_nt);

       BWExactSearchVectorForward(code_seq, start, end, 0, index->h_Oi.siz - 2,
				  ki1, li1, &index->h_C, &index->h_C1, &index->h_Oi,
				  &last_ki1, &last_li1, &forw1_nt);
     
       r_list.num_results = 0;
       r_list.read_index = 0;

       BWSearch1(code_seq, start, end, k1, l1, ki1, li1, 
		 &index->h_C, &index->h_C1, &index->h_O, &index->h_Oi, &r_list);

       for (size_t ii = 0; ii < r_list.num_results; ii++) {
	 r = &r_list.list[ii];

	 if(!r->num_mismatches)
	   error = 0;
	 else
	   error = r->err_kind[0];

	 pos = r->err_pos[0];  

	 direction = r->dir;

	 len_calc = len;
	 if (error == DELETION) {
	   len_calc--;
	 } else if (error == INSERTION) {
	   len_calc++;
	 }

	 // generating cigar
	 sprintf(quality_clipping, "%i", NONE_HARD_CLIPPING);
	 if (error == 0) {
	   sprintf(cigar, "%lu=\0", len);
	   num_cigar_ops = 1;
	   memcpy(seq_dup, seq_strand, len);
	   seq_dup[len] = '\0';
		 
	 } else if (error == MISMATCH) {
	   if (pos == 0) {
	     //Positive strand
	     if(type) { 
	       sprintf(cigar, "1S%luM\0", len-1); 
	       start_mapping++;
	     }
	     else { 
	       sprintf(cigar, "%luM1S\0", len-1); 
	     }
	     num_cigar_ops = 2;
	   } else if (pos == len - 1) {
	     //Positive strand
	     if(type) { 
	       sprintf(cigar, "%luM1S\0", len - 1); 
	     }
	     else{ 
	       sprintf(cigar, "1S%luM\0", len-1); 
	       start_mapping++;
	     }
	     num_cigar_ops = 2;
	   } else {
	     sprintf(cigar, "%luM\0", len);
	     num_cigar_ops = 1;
	   }
	   memcpy(seq_dup, seq_strand, len);
	   seq_dup[len] = '\0';
	   //printf("MISMATCH\n");
	 } else if (error == INSERTION) {
	   //printf("INSERTION\n");
	   if (pos == 0) {
	     if(type) {
	       sprintf(cigar, "1M1D%luM\0", len - 1); 
	     }
	     else{ 
	       sprintf(cigar, "%luM1D1M\0", len - 1); 
	     }	      
	   } else if (pos == len - 1) {
	     if(type) { 
	       sprintf(cigar, "%luM1D1M\0", len - 1); 
	     }
	     else{ 
	       sprintf(cigar, "1M1D%luM\0", len - 1); 
	     }
	   } else {
	     if(type) {
	       if(r->dir)
		 sprintf(cigar, "%iM1D%luM\0", pos, len - pos);
	       else
		 sprintf(cigar, "%iM1D%luM\0", pos + 1, len - pos - 1);
	     } else { 
	       if(r->dir)
		 sprintf(cigar, "%luM1D%dM\0", len - pos, pos);
	       else
		 sprintf(cigar, "%luM1D%dM\0", len - pos - 1, pos + 1);
	     }
	   }
	   num_cigar_ops = 3;
	   memcpy(seq_dup, seq_strand, len);
	   seq_dup[len] = '\0';
	 } else if (error == DELETION) {	     
	   //printf("DELETION\n");
	   if (pos == 0) {
	     if(type) { 
	       sprintf(cigar, "1I%luM\0", len -1); 
	     } else{ 
	       sprintf(cigar, "%luM1I\0", len -1); 
	       start_mapping++;
	     }
	     num_cigar_ops = 2;		
	   } else if (pos == len - 1) {
	     if(type) { 
	       sprintf(cigar, "%luM1I\0", len -1); 
	       start_mapping++;
	     } else{ 
	       sprintf(cigar, "1I%luM\0", len -1); 
	     }
	     num_cigar_ops = 2;
	   } else {
	     if(type) { 
	       sprintf(cigar, "%dM1I%luM\0", pos, len - pos - 1); 
	     }
	     else{ 
	       sprintf(cigar, "%luM1I%dM\0", len - pos - 1, pos); 
	     }
	     num_cigar_ops = 3;
	   }
	   memcpy(seq_dup, seq_strand , len );
	   seq_dup[len] = '\0';
		 
	 }else{
	   //printf("NUM MAPPINGS %lu -> POS %d -> ERROR %d -> (%lu):%s", num_mappings, pos, error, len, seq);
	   continue;
	 }
	 k_start = r->k;
	 l_start = r->l;

	 tot_alignments += (l_start - k_start);
	 // check filter_read_mappings for bisulfite case
	 if (tot_alignments >  bwt_optarg->filter_read_mappings) {
	   filter_exceeded = 1;
	   LOG_DEBUG_F("Filter exceeded: num. read mappings: %i (total: %i > filter %i)\n", 
		       l_start - k_start, tot_alignments, bwt_optarg->filter_read_mappings);
	   break;
	 }

	 for (size_t j = k_start; j <= l_start; j++) {
	   if (index->S.ratio == 1) {
	     key = (direction)
	       ? index->Si.siz - index->Si.vector[j] - len_calc - 1
	       : index->S.vector[j];
	   } else {
	     key = (direction)
	       ? index->Si.siz - getScompValue(j, &index->Si, &index->h_C,
					       &index->h_Oi) - len_calc - 1
	       : getScompValue(j, &index->S, &index->h_C, &index->h_O);
	   }
	   idx = binsearch(index->karyotype.offset, index->karyotype.size, key);
	   if(key + len_calc <= index->karyotype.offset[idx]) {
	     start_mapping = index->karyotype.start[idx-1] + (key - index->karyotype.offset[idx-1]);
	     // save all into one alignment structure and insert to the list
	     alignment = alignment_new();

	     alignment_init_single_end(NULL, strdup(seq_dup), strdup(quality_clipping), type, 
				       idx - 1,
				       start_mapping,
				       strdup(cigar), num_cigar_ops, 254, 1, (num_mappings > 0), 0, NULL, alignment);
	     //printf("\n\nStart = %lu\nEnd   = %lu\nDif   = %lu\n", start, end, end - start);
	     //alignment_print(alignment);
	     array_list_insert((void*) alignment, tmp_mapping_list);
	   }
	 }//end for k and l
       }//end for

       if (filter_exceeded) {
	 //printf("Limit Exceeded %d\n", bwt_optarg->filter_read_mappings);
	 array_list_clear(tmp_mapping_list, (void *)alignment_free);
	 if (!mapping_list->size) {
	   array_list_set_flag(BS_EXCEEDED, mapping_list);
	 }
	 goto exit;
       }
     
       //Search for equal BWT mappings and set the mappings that will be delete
       int n_mappings = array_list_size(tmp_mapping_list);
       alignment_t *alig_1, *alig_2;
       unsigned int *delete_mark = (unsigned int *)calloc(n_mappings, sizeof(unsigned int));
       const int max_distance = 10;
     
       for (int a1 = n_mappings - 1; a1 >= 1; a1--) {
	 if (!delete_mark[a1]) {
	   alig_1 = array_list_get(a1, tmp_mapping_list);
	   for (int a2 = a1 - 1; a2 >= 0; a2--) {
	     alig_2 = array_list_get(a2, tmp_mapping_list);
	     size_t dist = abs(alig_1->position - alig_2->position);
	     if (alig_1->chromosome == alig_2->chromosome &&
		 dist < max_distance && 
		 !delete_mark[a2]) {
	       //Same chromosome && same position
	       if (alig_1->num_cigar_operations < alig_2->num_cigar_operations) {
		 delete_mark[a2] = 1;
	       } else {
		 delete_mark[a1] = 1;
	       }
	     }
	   }
	 }
       }
     
       if (array_list_get_flag(mapping_list) == 1 &&
	   n_mappings) {
	 array_list_clear(mapping_list, (void *)bwt_anchor_free);
       }

       //Delete all set mappings
       int primary_delete = 0;
       int header_len;
       for (int m = n_mappings - 1; m >= 0; m--) {
	 alig_1 = array_list_remove_at(m, tmp_mapping_list);
	 if (delete_mark[m]) {
	   if (!is_secondary_alignment(alig_1)) { primary_delete = 1; }
	   alignment_free(alig_1);
	 } else {
	   set_secondary_alignment(num_mappings > 0, alig_1);
	   num_mappings++;
	   header_len = strlen(read->id);
	   alig_1->query_name = (char *) malloc(sizeof(char) * (header_len + 1));
	   get_to_first_blank(read->id, header_len, alig_1->query_name);
	   bwt_cigar_cpy(alig_1, read->quality);
	   alig_1 = add_optional_fields(alig_1, n_mappings);
	   array_list_insert((void*) alig_1, mapping_list);
	 }
       }

       if (primary_delete) { 
	 alig_1 = array_list_get(0, mapping_list); 
	 set_secondary_alignment(0, alig_1);
       }
       free(delete_mark);

       if (n_mappings == 0) { // no alignments found
	 //printf("Not mapping, maybe anchors (flag %i [%lu])\n", array_list_get_flag(mapping_list_tmp), array_list_size(mapping_list_tmp));
	 //printf("Size %lu\n", array_list_size(mapping_list_tmp));
	 if (array_list_get_flag(mapping_list_tmp) == BS_ANCHORS) {
	   //printf("Not mapping, anchors before [%lu]\n", array_list_size(mapping_list_tmp));
	   array_list_t *forward_anchor_list, *backward_anchor_list;
	   int new_type = !type; //1 == (+)

	   //printf("BACKWARD (+)\n");
	   //printf("fowr1 %i : back1 %i\n", forw1_nt, back1_nt);
	   size_t num_elem;
	   if (back1_nt > min_size) {
	     LOG_DEBUG("1 INIT ANCHOR LIST");
	     __bwt_generate_anchor_list(last_k1, last_l1, back1_nt, bwt_optarg, 
					index, new_type, anchor_list, BACKWARD_ANCHOR, 0);
	     // convert anchors in anchor_list, into cals in mapping_list_tmp
	     LOG_DEBUG("1 END ANCHOR LIST");
	     num_elem = array_list_size(anchor_list);
	     for (int elem = num_elem - 1; elem >= 0; elem--) {
	       anchor = array_list_remove_at(elem, anchor_list);
	       cal = convert_bwt_anchor_to_CAL(anchor, end - back1_nt - 1, end - 1);
	       array_list_insert(cal, mapping_list_tmp);
	     }
	     // necessary??
	     array_list_clear(anchor_list, (void *) bwt_anchor_free);
	   } else {
	     back1_nt = 0;
	   }
	   //printf("FORWARD (+)\n");
	   if (forw1_nt > min_size) {
	     LOG_DEBUG("2 INIT ANCHOR LIST");
	     __bwt_generate_anchor_list(last_ki1, last_li1, forw1_nt, bwt_optarg, 
					index, new_type,  mapping_list_tmp, FORWARD_ANCHOR, 0);
	     LOG_DEBUG("2 END ANCHOR LIST");
	     // convert anchors in anchor_list, into cals in mapping_list_tmp
	     num_elem = array_list_size(anchor_list);
	     for (int elem = num_elem - 1; elem >= 0; elem--) {
	       anchor = array_list_remove_at(elem, anchor_list);
	       cal = convert_bwt_anchor_to_CAL(anchor, start, start + forw1_nt - 1);
	       array_list_insert(cal, mapping_list_tmp);
	     }
	     // necessary??
	     array_list_clear(anchor_list, (void *) bwt_anchor_free);
	   } else {
	     forw1_nt = 0;
	   }
	   //printf("Size %lu\n", array_list_size(mapping_list_tmp));

	   // jump
	   start+= (jump_error + forw1_nt);
	   end  -= (jump_error + back1_nt);
	 } else {
	   // there is no alignments nor anchors
	   // skip some positions, and try again
	   start+=jump_error;
	   end-=jump_error;
	 }
	 //printf("Size %lu\n", array_list_size(mapping_list_tmp));
       } else { // alignments found

	 // insert the new mappings in the output array
	 n_mappings = array_list_size(mapping_list);
	 if (first_loop == 0) {
	   //printf("Mapping in the first loop [%lu]\n", array_list_size(mapping_list_tmp));
	   // si la lista destino tiene anchors, se eliminan
	   if (array_list_get_flag(mapping_list_tmp) == BS_ANCHORS) {
	     //printf("Mappings, clear anchors\n");
	     array_list_clear(mapping_list_tmp, (void *) bwt_anchor_free);
	   }
	   array_list_set_flag(BS_ALIGNMENTS, mapping_list_tmp);
	   for (int m = n_mappings - 1; m >= 0; m--) {
	     alig_1 = array_list_remove_at(m, mapping_list);
	     array_list_insert((void*) alig_1, mapping_list_tmp);
	   }
	 } else {
	   //printf("Mappings not the first loop\n");
	   array_list_set_flag(BS_ANCHORS, mapping_list_tmp);
	   // convertir alignments en bwt_anchor_t
	   for (int m = n_mappings - 1; m >= 0; m--) {
	     alig_1 = array_list_remove_at(m, mapping_list);
	     bwt_anchor_t *anchor_tmp = bwt_anchor_new(alig_1->seq_strand,
						       alig_1->chromosome,
						       alig_1->position + start,
						       alig_1->position + end - start,
						       0);
	     array_list_insert((void*) anchor_tmp, mapping_list_tmp);
	     alignment_free(alig_1);
	   }
	 }
	 // found alignment
	 start = end;
       }
       //printf("\n");

       first_loop++;
       array_list_clear(mapping_list, (void *) NULL);
       array_list_clear(tmp_mapping_list, (void *) NULL);
       ////////////////////////////////
     }//end while
     ////////////////////////////////

 exit:
     free(r_list.list);
     free(code_seq);
     free(seq_strand);
     free(k1);
     free(l1);
     free(ki1);
     free(li1);

     free(seq_dup);
     free(quality_clipping);
     array_list_free(tmp_mapping_list, (void *)NULL);
     array_list_free(mapping_list, (void *)NULL);
     array_list_free(anchor_list, (void *)bwt_anchor_free);

     return array_list_size(mapping_list_tmp);
}

//-----------------------------------------------------------------------------
