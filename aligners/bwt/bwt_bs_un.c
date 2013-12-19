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
    array_list_set_flag(BS_ALIGNMENTS, mapping_list_tmp);
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

  array_list_t *region_list = array_list_new(bwt_optarg->filter_read_mappings, 1.25f,
					     COLLECTION_MODE_ASYNCHRONIZED);

  //cal_t *cal;
  region_t *region;
  bwt_anchor_t * anchor;
     
  int min_size = 20;
  size_t first_loop = 0;
  size_t id = 0;
  ///////////////////////////////
  while ((int)end - (int)start >= min_size) {
    ///////////////////////////////
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
      if (!mapping_list_tmp->size) {
	size_t flag = array_list_get_flag(mapping_list_tmp);
	if (flag == BS_ALIGNMENTS)
	  array_list_clear(mapping_list_tmp, (void *)alignment_free);
	else if (flag == BS_ANCHORS)
	  array_list_clear(mapping_list_tmp, (void *)region_bwt_free);
	else
	  array_list_clear(mapping_list_tmp, (void *)NULL);

	array_list_set_flag(BS_EXCEEDED, mapping_list_tmp);
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
     
    if (array_list_get_flag(mapping_list) == BS_ANCHORS && n_mappings) {
      array_list_clear(mapping_list_tmp, (void *) region_bwt_free);
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
      if (array_list_get_flag(mapping_list_tmp) == BS_ANCHORS) {
	int new_type = !type; //1 == (+)

	size_t num_elem;
	if (back1_nt > min_size) {
	  LOG_DEBUG("1 INIT ANCHOR LIST");
	  __bwt_generate_anchor_list(last_k1, last_l1, back1_nt, bwt_optarg, 
				     index, new_type, anchor_list, BACKWARD_ANCHOR, 0);
	  // convert anchors in anchor_list, into cals in mapping_list_tmp
	  LOG_DEBUG("1 CONVERT ANCHOR LIST");
	  num_elem = array_list_size(anchor_list);
	  for (int elem = num_elem - 1; elem >= 0; elem--) {
	    anchor = array_list_remove_at(elem, anchor_list);
	    region = region_bwt_new(anchor->chromosome, anchor->strand,
				    anchor->start, anchor->end,
				    end - back1_nt - 1, end - 1,
				    read->length, id++);
	    array_list_insert(region, mapping_list_tmp);
	  }
	  // necessary??
	  //array_list_clear(anchor_list, (void *) bwt_anchor_free);
	} else {
	  back1_nt = 0;
	}
	//printf("FORWARD (+)\n");
	if (forw1_nt > min_size) {
	  LOG_DEBUG("2 INIT ANCHOR LIST");
	  __bwt_generate_anchor_list(last_ki1, last_li1, forw1_nt, bwt_optarg, 
				     index, new_type,  anchor_list, FORWARD_ANCHOR, 0);
	  LOG_DEBUG("2 CONVERT ANCHOR LIST");
	  // convert anchors in anchor_list, into cals in mapping_list_tmp
	  num_elem = array_list_size(anchor_list);
	  for (int elem = num_elem - 1; elem >= 0; elem--) {
	    anchor = array_list_remove_at(elem, anchor_list);
	    region = region_bwt_new(anchor->chromosome, anchor->strand,
				    anchor->start, anchor->end,
				    start, start + forw1_nt - 1,
				    read->length, id++);
	    array_list_insert(region, mapping_list_tmp);
	  }
	  // necessary??
	  //array_list_clear(anchor_list, (void *) bwt_anchor_free);
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
	  array_list_clear(mapping_list_tmp, (void *) region_bwt_free);
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
	  region = region_bwt_new(alig_1->chromosome, alig_1->seq_strand,
				  alig_1->position, alig_1->position + end - start,
				  start, start + forw1_nt - 1,
				  read->length, id++);
	  array_list_insert((void*) region, mapping_list_tmp);
	  alignment_free(alig_1);
	}
      }
      // found alignment
      start = end;
    }
    //printf("\n");

    first_loop++;
    array_list_clear(region_list, (void *) region_bwt_free);
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
  array_list_free(region_list, (void *)region_bwt_free);

  return array_list_size(mapping_list_tmp);
}

//-----------------------------------------------------------------------------

size_t bwt_generate_cals_bs_un(char *seq, char *seq2, size_t seed_size, bwt_optarg_t *bwt_optarg, 
			       bwt_index_t *index, bwt_index_t *index2, array_list_t *cal_list) {

  array_list_t *mapping_list = array_list_new(100, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  const int nstrands = 2;
  const int nchromosomes = 30; //TODO: Parameter
  
  linked_list_t ***cals_list = (linked_list_t ***)malloc(sizeof(linked_list_t **)*nstrands);

  for (unsigned int i = 0; i < nstrands; i++) {
    cals_list[i] = (linked_list_t **)malloc(sizeof(linked_list_t *)*nchromosomes);
    for (unsigned int j = 0; j < nchromosomes; j++) {
      cals_list[i][j] = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    }
  }

  int seed_id = 0;
  size_t len = strlen(seq);
  size_t len2 = strlen(seq2);

  size_t offset, num_seeds;
  
  size_t num_mappings;
  size_t start, end, offset_end = len - 16;
  size_t n_seeds = num_seeds;
  size_t offset_inc = ceil(1.0f * len / (num_seeds + 1));
  if (offset_inc <= 0) offset_inc = 1;

  /*
  num_mappings = bwt_map_exact_seed_bs(code_seq, len, start, end - 1,
				       bwt_optarg, index,  mapping_list, seed_id++, 0);
  insert_seeds_and_merge(mapping_list, cals_list, len);
  */
  LOG_DEBUG("Insert seed and merge");
  insert_seeds_and_merge(cal_list, cals_list, len);
  array_list_clear(cal_list, (void *) region_bwt_free);
  //Store CALs in Array List for return results
  size_t start_cal, end_cal;
  size_t seq_start, seq_end;
  seed_region_t *s, *s_first, *s_last, *s_aux, *seed_region;
  cal_t *cal;
  short_cal_t *short_cal, *short_cal_aux;
  linked_list_iterator_t itr, itr2;
  linked_list_item_t *list_item_cal, *list_item_aux;
  size_t min_cal_size = 20;//cal_optarg->min_cal_size; TODO:PARAMETER
  const int max_intron_size = 500000; //TODO: Parameter
  size_t read_length = len;

  for (unsigned int j = 0; j < nchromosomes; j++) {
    for (unsigned int i = 0; i < nstrands; i++) {
      linked_list_iterator_init(cals_list[i][j], &itr);
      list_item_cal = linked_list_iterator_list_item_curr(&itr);
      while ((list_item_cal != NULL )) {
	short_cal = (short_cal_t *)list_item_cal->item;
	if (short_cal->end - short_cal->start + 1 >= min_cal_size) {
	  linked_list_t *list_aux = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	  while (s = (seed_region_t *)linked_list_remove_last(short_cal->sr_list)) {
	    //TODO: Change all parameters to seed_region_t
	    append_seed_region_linked_list(list_aux,
					   s->read_start, s->read_end, 
					   s->genome_start, s->genome_end, 
					   s->id);	    
	    seed_region_free(s);	    
	  }
	  cal = cal_new(j, i, short_cal->start, 
			short_cal->end, short_cal->num_seeds, 
			list_aux, short_cal->sr_duplicate_list);
 	  array_list_insert(cal, cal_list);
	  
	  short_cal->sr_duplicate_list = NULL;
	  
	  //============= Search if one seed map near to this CAL ================//
	  s_first = (seed_region_t *)linked_list_get_first(list_aux);	  
	  s_last = (seed_region_t *)linked_list_get_last(list_aux);
	  
	  if (s_first && s_first->read_start > seed_size) {
	    //Search start seed <-----
	    list_item_aux = list_item_cal->prev;
	    while (list_item_aux) {
	      short_cal_aux = list_item_aux->item;
	      if (!short_cal_aux ||
		  short_cal_aux->end - short_cal_aux->start + 1 >= min_cal_size || 
		  short_cal->start - short_cal_aux->end >= max_intron_size) {
		break;
	      }
	   
	      linked_list_iterator_init(short_cal_aux->sr_list, &itr2);
	      s = linked_list_iterator_curr(&itr2);
	      while (s != NULL) {
		if (s->read_end < s_first->read_start) {
		  seed_region = seed_region_new(s->read_start, s->read_end,
						s->genome_start, s->genome_end, 0);	
		  array_list_insert(seed_region, cal->candidates_seeds_start);
		}
		s = linked_list_iterator_next(&itr2);

	      }
	      list_item_aux = list_item_aux->prev;	      
	    }
	  }	

	  if (s_last && s_last->read_end < read_length - seed_size) {
	    // Search start seed ----->
	    list_item_aux = list_item_cal->next;
	    while (list_item_aux) {
	      short_cal_aux = list_item_aux->item;
	      if (!short_cal_aux ||
		  short_cal_aux->end - short_cal_aux->start + 1 >= min_cal_size || 
		  short_cal_aux->start - short_cal->end >= max_intron_size) {
		break;
	      }

	      linked_list_iterator_init(short_cal_aux->sr_list, &itr2);
	      s = linked_list_iterator_curr(&itr2);
	      while (s != NULL) {
		if (s->read_start > s_last->read_end) {
		  seed_region = seed_region_new(s->read_start, s->read_end,
						s->genome_start, s->genome_end, 1);	
		  array_list_insert(seed_region, cal->candidates_seeds_end);
		}
		s = linked_list_iterator_next(&itr2);
	      }
	      list_item_aux = list_item_aux->next;
	    }	    
	  }
	  //======================================================================//
        } 
	linked_list_iterator_next(&itr);
	list_item_cal = linked_list_iterator_list_item_curr(&itr);
      }
    }
  }

  for (unsigned int i = 0; i < nstrands; i++) {
    for (unsigned int j = 0; j < nchromosomes; j++) {
      linked_list_free(cals_list[i][j], (void *)short_cal_free);
    }
    free(cals_list[i]);
  }
  
  array_list_free(mapping_list, NULL);
  free(cals_list);

  return array_list_size(cal_list);
}

//-----------------------------------------------------------------------------
