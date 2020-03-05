/* Function for finding edge loops, duplicates and bi-directional edges.
 *
 * Author: Alex T. Kalinka (alex.t.kalinka@gmail.com)
 *
 * Will be loaded as a shared object into R (dyn.load('name.so')) and called from within R.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <vector>
#include <iterator>
#include <algorithm>
#include <stdint.h>

#include <R.h>

extern "C" {

using namespace std;


void edgeDuplicates(int *edgeA, int *edgeB, int *numedg, int *loops, int *dups, bool *verbose) 

	{

	int i, sum;
	float prog;
	uint_fast32_t tpair; // int_fast32_t uses 32-bit at least, but will use 64-bit if possible.
	uint_fast32_t run_max = 0;
	vector<uint_fast32_t> pair_ids; // Unique IDs for each edge using Cantor's pairing function (see: http://en.wikipedia.org/wiki/Pairing_function ).
	vector<uint_fast32_t>::iterator pit;

	// Loop through edges and ask if they're loops, duplicates, or bidirectional.
	for(i = 0; i < *numedg; i++){
		
		if(*verbose){
			prog = (i+0.0)/(*numedg-1)*100;

			Rprintf("\r   Checking for loops and duplicate edges... %3.3f%%",prog);

			R_FlushConsole();
			R_ProcessEvents();
			}

		if(edgeA[i] == edgeB[i]){ // Loop.
			loops[i] = 1;
			continue;
		}else{
			// Cantor's pairing function.
			sum = edgeA[i] + edgeB[i];
			if(edgeA[i] < edgeB[i]){
				tpair = (sum*(sum+1)/2) + edgeA[i];
			}else{
				tpair = (sum*(sum+1)/2) + edgeB[i];
				}

			if(tpair <= run_max){
				pit = find(pair_ids.begin(), pair_ids.end(), tpair );
			}else{
				pit = pair_ids.end();
				run_max = tpair;
				}
			if( pit != pair_ids.end() ){ // Duplicate or bi-directional edge.
				dups[i] = 1;
			}else{
				pair_ids.push_back( tpair );
				}

			}

		}

	}


	}





