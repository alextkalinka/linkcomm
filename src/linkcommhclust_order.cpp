/* Code to order hierarchical clusters horizontally for plotting in R as an object of class "hclust".
 *
 * Author: Alex T. Kalinka (alex.t.kalinka@gmail.com)
 *
 * Will be loaded as a shared object into R (dyn.load('name.so')) and called from within R using .C("hclustPlotOrder",...). */

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>

#include <R.h>

extern "C" {

using namespace std;

void hclustPlotOrder(int *numedg, int *hca, int *hcb, int *order)
	{
	
	vector<int> rowM;
	vector<int> colM;
	vector<int> done;
	vector<int> old;
	int sub = 0, num = *numedg-1, next = -1;

	copy(hca, hca + (*numedg-1), back_inserter(rowM));
	copy(hcb, hcb + (*numedg-1), back_inserter(colM));

	/* Start with highest cluster and work down.
	 * Tightest clusters positioned on the left.
	 * The tightest cluster is a singleton. */ 
	while(num > 0){
		if(find(done.begin(),done.end(),num-1) != done.end()){
			num--;
			continue;
			}
		if(rowM.at(num-1) < 0){
			order[sub] = abs(rowM.at(num-1));
			sub++;
			if(colM.at(num-1) < 0){
				order[sub] = abs(colM.at(num-1));
				sub++;
			}else{	
				next = colM.at(num-1) - 1;
				}
		}else{
			old.push_back( colM.at(num-1) - 1 );
			next = rowM.at(num-1) - 1;
			}
		
		while(next >= 0){
			done.push_back(next);
			if(rowM.at(next) < 0){
				order[sub] = abs(rowM.at(next));
				sub++;
				if(colM.at(next) < 0){
					order[sub] = abs(colM.at(next));
					sub++;
					if(old.size() > 0){
						next = old.back();
						old.pop_back();
					}else{
						next = -1;
						}
				}else{
					next = colM.at(next) - 1;
					}
			}else{
				old.push_back( colM.at(next) - 1 );
				next = rowM.at(next) - 1;
				}
			}
		num--;
		}

	}

	}
			




