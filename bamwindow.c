// bamwindow.c - Print the number of reads aligning to windows of a fixed size in a BAM file.
// Copyright (C) 2015 Alastair Droop, The Leeds MRC Medical Bioinformatics Centre <a.p.droop@leeds.ac.uk>

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include "sam.h"

#define PROG_NAME "bamwindow"
#define PROG_VERSION "2.1a (2019-02-22)"

// Global variables:
char print_empty = 0;
char match_mode = 0;
char *target_name = "";
unsigned long window_start = 0;
unsigned long window_end = 0;
unsigned long window_count = 0;

// -m0, also the default:
static int pileup_overlap(const bam1_t *b, void *data){
  window_count ++;
  return 0;
}

// -m1:
static int pileup_start(const bam1_t *b, void *data){
  uint32_t readstart = b->core.pos;
  if(readstart < window_start) return 0;
  if(readstart >= window_end) return 0;
  window_count ++;
  return 0;
}

// -m2:
static int pileup_midpoint(const bam1_t *b, void *data){
    uint32_t midpoint = b->core.pos + ((bam_endpos(b) - b->core.pos) / 2);
    if(midpoint < window_start) return 0;
    if(midpoint >= window_end) return 0;
    window_count ++;
    return 0;
}

typedef void (*pileup_callback)(const bam1_t *b, void *data);
pileup_callback callback;

void print_usage(){
    fprintf(stdout, "usage: %s [-hvem] n file [region]\n", PROG_NAME);
}

void print_version(){
    fprintf(stdout, "%s %s (HTSlib version %s)\n", PROG_NAME, PROG_VERSION, hts_version());
}

void print_help(){
    printf("Print the number of reads aligning to windows of a fixed size\n");
    printf("in a BAM file.\n\n");
    print_usage();
    printf("\nIf no region is specified, for each target (usually a chromosome)\n");
    printf("in the BAM file, windows are created starting at the beginning of the\n");
    printf("target and incrementing by n. The last window is clipped to to the end\n");
    printf("of the target. If a region is specified, windows are created starting\n");
    printf("at the start of the region, and incrementing by n. The last window is\n");
    printf("clipped to the end of the region or the end of the target. Regions are\n");
    printf("1-based, so both the start and end coordinate are included in a region.\n");
    printf("The read start is left-most read nucleotide.\n");
    printf("For midpoint mapping, the read length is calculated excluding any clipped\n");
    printf("regions. Return ranges are 1-based and inclusive.\n");
    printf("\nThe match mode defines how reads are matched to windows:\n");
    printf(" -m0: match by overlap (reads can therefore be counted multiple times)\n");
    printf(" -m1: match by read start\n");
    printf(" -m2: match by read midpoint\n");
    printf("\nOptions & arguments:\n");
    printf("-h     : print this help and quit\n");
    printf("-v     : print the version and quit\n");
    printf("-e     : print empty regions\n");
    printf("-m     : the matching mode [0, 1 or 2] (default is 0)\n");
    printf("n      : the window size\n");
    printf("file   : BAM file to process\n");
    printf("region : A region (target:start-end)\n");
    printf("\nOutput: <target> <start> <end> <count>.\n");
}

void process_region(samFile *in, bam1_t *b, hts_idx_t *index, bam_hdr_t *header, char *range){
    if (window_start - window_end < 1) return;
    int result = 0;
    hts_itr_t *range_iter = NULL;
    window_count = 0;
    range_iter = sam_itr_querys(index, header, range);
    while ((result = sam_itr_next(in, range_iter, b)) >= 0) callback(b, NULL);
    if (window_count == 0 && print_empty == 0) return;
    fprintf(stdout, "%s\t%lu\t%lu\t%lu\n", target_name, window_start + 1, window_end, window_count);
}

int main(int argc, char *argv[]){
    int option = 0;
    int exit_code = 0;
    uint32_t window_size;
    int32_t target_number;
    uint32_t target_end = 0;
    samFile *input_file = 0;
    bam_hdr_t *input_header = NULL;
    bam1_t *bam_data = NULL;
    hts_idx_t *input_index = NULL;
    uint32_t window_length = 0;
    int current_max_rlen = 100;
    char *range_string = NULL;
    
    //Process options:
    while((option = getopt(argc, argv, "hvem:")) >= 0){
        switch(option){
            case 'h':
                print_help();
                return 0;
            case 'v':
                print_version();
                return 0;
            case 'e': print_empty = 1; break;
            case 'm':
              if(*optarg == '0') match_mode = 0;
              else if(*optarg == '1') match_mode = 1;
              else if(*optarg == '2') match_mode = 2;
              else{
                printf("match mode must be 0, 1 or 2\n");
                return 1;
              }
              break;
        }
    }
    
    //Allocate target string:
    range_string = (char*)malloc(current_max_rlen * sizeof(char));
    if (range_string == NULL){
        fprintf(stderr, "ERROR: Failed to reallocate target string\n");
        exit_code = 1;
        goto cleanup;
    }
    
    // Set the pileup callback:
    if(match_mode == 0) callback = (void *)pileup_overlap;
    else if (match_mode == 1) callback = (void *)pileup_start;
    else if (match_mode == 2) callback = (void *)pileup_midpoint;

    //Make sure we have 2 or 3 arguments left:
    if ((argc - optind) < 2  || (argc - optind) > 3 ){
        print_usage();
        return 1;
    }

    //Extract the window size:
    window_size = atoi(argv[optind]);
    optind ++;
    if(window_size < 1){
      printf("invalid window size\n");
      return 1;
    }

    //Open the file specified:
    if ((input_file = sam_open(argv[optind], "r")) == 0) {
        printf("ERROR: Failed to open file \"%s\"\n", argv[optind]);
        exit_code = 1;
        goto cleanup;
    }

    //Read in the file header:
    if ((input_header = sam_hdr_read(input_file)) == 0) {
        fprintf(stderr, "[main_samview] fail to read the header from \"%s\".\n", argv[optind]);
        exit_code = 1;
        goto cleanup;
    }

    // Read in the file index:
    input_index = sam_index_load(input_file, argv[optind]); // load index
    if (input_index == 0) { // index is unavailable
        fprintf(stderr, "ERROR: Failed to open file index\n");
        exit_code = 1;
        goto cleanup;
    }
    
    // Initialise the BAM data object:
    bam_data = bam_init1();
    
    optind++;
    if (argc > optind){
        // Process only the given range:
        hts_itr_t *iter = sam_itr_querys(input_index, input_header, argv[optind]);
        target_name = input_header->target_name[iter->tid];
        target_end = iter->end;
        window_start = iter->beg;
        int target_max_rlen = strlen(target_name) + (ceil(log10(target_end)) * 2) + 3;
        if (target_max_rlen > current_max_rlen){
            current_max_rlen = target_max_rlen;
            range_string = (char*)realloc(range_string, current_max_rlen * sizeof(char));
            if (range_string == NULL){
                fprintf(stderr, "ERROR: Failed to reallocate target string\n");
                exit_code = 1;
                goto cleanup;
            }
        }
        while(1){
            window_end = window_start + window_size;
            if (window_end > target_end) window_end = target_end;
            window_length = window_end - window_start;
            sprintf(range_string, "%s:%lu-%lu", target_name, window_start, window_end);
            process_region(input_file, bam_data, input_index, input_header, range_string);
            if (window_length < window_size) break;
            window_start = window_end;
        }
    } else {
        //Run through the whole file:
        for (target_number = 0; target_number < input_header->n_targets; target_number++){
            target_name = input_header->target_name[target_number];
            target_end = input_header->target_len[target_number];
            window_start = 0;
            int target_max_rlen = strlen(target_name) + (ceil(log10(target_end)) * 2) + 3;
            if (target_max_rlen > current_max_rlen){
                current_max_rlen = target_max_rlen;
                range_string = (char*)realloc(range_string, current_max_rlen * sizeof(char));
                if (range_string == NULL){
                    fprintf(stderr, "ERROR: Failed to reallocate target string\n");
                    exit_code = 1;
                    goto cleanup;
                }
            }
            while (1){
                window_end = window_start + window_size;
                if (window_end > target_end) window_end = target_end;
                window_length = window_end - window_start;
                sprintf(range_string, "%s:%lu-%lu", target_name, window_start, window_end);
                process_region(input_file, bam_data, input_index, input_header, range_string);
                if (window_length < window_size) break;
                window_start = window_end;
            }
        }
    }

    //Cleanup:
cleanup:
    bam_destroy1(bam_data);
    bam_hdr_destroy(input_header);
    sam_close(input_file);
    hts_idx_destroy(input_index);
    if (range_string != NULL) free(range_string);
    return exit_code;
}
