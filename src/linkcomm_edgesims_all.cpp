/* Code to calculate edge similarities for hierarchical clustering of edges to produce link communities
 *  when including all edges, not just those that share a node.
 * For undirected networks Jaccard coefficients will be calculated, for directed and/or weighted networks the Tanimoto coefficient.
 *
 * Author: Alex T. Kalinka (alex.t.kalinka@gmail.com)
 *
 * Will be loaded as a shared object into R (dyn.load('name.so')) and called from within R.
 * Based on the Link communities derived from the algorithm in:
 * Ahn et al. (2010). Link communities reveal multiscale complexity in networks. Nature 466:761-765. 
 *

*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <set>
#include <map>
#include <numeric>

#include <R.h>

extern "C" {

using namespace std;


void nodeNeighbourhood(map<int, set<int> > &nN, map<int, map<int,double> > &wM, vector<int> &eA, vector<int> &eB, int node, unsigned int *numedg, bool *weighted, double *weights)

	{

	unsigned int finda = 0, findb = 0, add;
	bool end = FALSE, first = TRUE;

	// Get first-order node neighbourhood, and separately their weights (if applicable).
	while(end == FALSE){
		if(finda != *numedg){
			if(first){add = 0;}else{add = finda+1;}
			if(add > *numedg){add = *numedg;}
			finda = find(eA.begin()+add, eA.end(), node) - eA.begin();
			if(finda != *numedg){
				nN[node].insert(eB.at(finda));
				if(*weighted){
					wM[node].insert( pair<int,double>(eB.at(finda), weights[finda]) );
					}
				}
			}
		if(findb != *numedg){
			if(first){add = 0;}else{add = findb+1;}
			if(add > *numedg){add = *numedg;}
			findb = find(eB.begin()+add, eB.end(), node) - eB.begin();
			if(findb != *numedg){
				nN[node].insert(eA.at(findb));
				if(*weighted){
					wM[node].insert( pair<int,double>(eA.at(findb), weights[findb]) );
					}
				}
			}

		if(finda == *numedg && findb == *numedg){
			end = TRUE;
		}else{
			first = FALSE;
			}

		}

	}



void getDirectedWeights_all(map<int,float> &dW, set<int> &comm, vector<int> &eA, vector<int> &eB, int nonSA, int nonSB, unsigned int *numedg, double *dirw)
	
	{

	unsigned int i, add = 0, findNA = 0, findNB = 0;
	bool matchA = FALSE, matchB = FALSE, B = FALSE, first = TRUE;
	vector<int> commV(comm.size());

	copy(comm.begin(),comm.end(),commV.begin());

	for(i = 0; i < comm.size(); i++){

		while(matchA == FALSE){
			if(first){add = 0;}else{add = findNA + 1;}
			if(add > *numedg){add = *numedg;}
			findNA = find(eA.begin() + add, eA.end(), nonSA) - eA.begin();
			if(findNA != *numedg){
				if(eB.at(findNA) == commV.at(i)){
					B = TRUE;
					matchA = TRUE;
					}
			}else{
				matchA = TRUE;
				}
			first = FALSE;
			}

		first = TRUE;
			
		while(matchB == FALSE){
			if(first){add = 0;}else{add = findNB + 1;}
			if(add > *numedg){add = *numedg;}
			findNB = find(eA.begin() + add, eA.end(), nonSB) - eA.begin();
			if(findNB != *numedg){
				if(eB.at(findNB) == commV.at(i)){
					if(B){
						dW.insert( pair <int,float>(commV.at(i), 1.0) );
					}else{
						dW.insert( pair <int,float>(commV.at(i), *dirw) );
						}
					matchB = TRUE;
					}
			}else{
				if(B){
					dW.insert( pair <int,float>(commV.at(i), *dirw) );
				}else{
					dW.insert( pair <int,float>(commV.at(i), 1.0) );
					}
				matchB = TRUE;
				}
			first = FALSE;
			}

		findNA = 0; findNB = 0;
		matchA = FALSE; matchB = FALSE; B = FALSE; first = TRUE;

		}

	}



void getEdgeSimilarities_all(int *ea, int *eb, unsigned int *numedg, unsigned int *numnodes, unsigned int *rowlen, double *weights, bool *directed, double *dirweight, bool *weighted, bool *disk, double *dissvec, bool *bipartite, bool *verbose)

	{

	unsigned int i, j, k = 0;
	int first = -1, last = -1, sum = 0, runn = 0;
	double dotprod, absA, absB, numerat, denom, distm;
	float prog;
	vector<double> row;
	vector<double> temp;
	vector<int> edgeA;
	vector<int> edgeB;
	vector<int> inds;
	vector<int> nonshared;
	set<int> nodesI;
	set<int> nodesJ;
	set<int> diff;
	set<int> neighbFocal;
	set<int> neighbOther;
	set<int> neighbA;
	set<int> neighbB;
	set<int> common;
	set<int> total;
	set<int>::iterator sit;
	map<int, set<int> > nodeNeighb;
	map<int, map<int,double> > weightMap;
	map<int,double> mapA;
	map<int,double> mapB;
	map<int,double> mapC;
	map<int,double> mapD;
	map<int,float> dirWeights;
	vector<double> aI;
	vector<double> aJ;

	copy(ea, ea + *numedg, back_inserter(edgeA));
	copy(eb, eb + *numedg, back_inserter(edgeB));

	ofstream outfile;

	if(*disk){
		remove("linkcomm_diss.txt");
		outfile.open("linkcomm_diss.txt", ios::out);
		if(! outfile.is_open()){
				Rprintf("\nERROR: can't open linkcomm_diss.txt!\n"); return;
				}

		outfile.precision(7);
		}

	// Get first-order node neighbourhood for all nodes and hold in memory.
	for(i = 1; i <= *numnodes; i++){
		
		if(*verbose){
			prog = (i+0.0)/(*numnodes)*100;

			Rprintf("\r   Calculating edge similarities for %d edges...1/2... %3.2f%%",*numedg,prog);

			R_FlushConsole();
			R_ProcessEvents();
			}

		nodeNeighbourhood(nodeNeighb, weightMap, edgeA, edgeB, i, numedg, weighted, weights);

		}


	// Loop through edges and calculate edge similarities.
	for(i = 0; i < *numedg-1; i++){

		if(*verbose){
			prog = (i+0.0)/(*numedg-2)*100;

			Rprintf("\r   Calculating edge similarities for %d edges...2/2... %3.2f%%",*numedg,prog);

			R_FlushConsole();
			R_ProcessEvents();
			}

		if(*disk){
			row.assign(*numedg-1-i, 1);
			}

		neighbA = nodeNeighb[edgeA.at(i)];
		neighbB = nodeNeighb[edgeB.at(i)];
		// The node neighbourhood for the focal edge.
		set_union(neighbA.begin(),neighbA.end(),neighbB.begin(),neighbB.end(), inserter(neighbFocal, neighbFocal.begin()));

		neighbA.clear();
		neighbB.clear();

		if(*weighted){
			mapA = weightMap[edgeA.at(i)];
			mapB = weightMap[edgeB.at(i)];
			if(!(*bipartite)){
				mapA.insert( pair<int,double>(edgeA.at(i), 1) );
				mapB.insert( pair<int,double>(edgeB.at(i), 1) );
				}
			}

		// Add focal nodes (inclusive neighbourhood set).
		if(!(*bipartite)){
			neighbFocal.insert(edgeA.at(i));
			neighbFocal.insert(edgeB.at(i));
			}

		// Loop through edges for which we have not yet calculated a similarity score and do so.
		for(j = (i+1); j <= *numedg-1; j++){

			neighbA = nodeNeighb[edgeA.at(j)];
			neighbB = nodeNeighb[edgeB.at(j)];
			// The node neighbourhood for the nodes in edge j.
			set_union(neighbA.begin(),neighbA.end(),neighbB.begin(),neighbB.end(), inserter(neighbOther, neighbOther.begin()));

			//Rprintf("\ninds.size %d\nedgeA.at(inds.at(j)) %d\nedgeB.at(inds.at(j)) %d\n",inds.size(),edgeA.at(inds.at(j)),edgeB.at(inds.at(j)));


			if(!(*bipartite)){ // Inclusive node neighbourhood only if not bipartite.
				neighbOther.insert(edgeA.at(j));
				neighbOther.insert(edgeB.at(j));
				}

			if(*weighted){
				mapC = weightMap[edgeA.at(j)];
				mapD = weightMap[edgeB.at(j)];
				if(!(*bipartite)){
					mapC.insert( pair<int,double>(edgeA.at(j), 1) );
					mapD.insert( pair<int,double>(edgeB.at(j), 1) );
					}
				}

			//for(sit = neighbA.begin(); sit != neighbA.end(); sit++){
			//	Rprintf("\nneighbA: %d\n",*sit);
			//	}
			//for(sit = neighbB.begin(); sit != neighbB.end(); sit++){
			//	Rprintf("\nneighbB: %d\n",*sit);
			//	}

			set_intersection(neighbFocal.begin(),neighbFocal.end(),neighbOther.begin(),neighbOther.end(), inserter(common, common.begin()));

			if(*weighted || *directed){
				set_union(neighbFocal.begin(),neighbFocal.end(),neighbOther.begin(),neighbOther.end(), inserter(total, total.begin()));
				}

			if(!*weighted && !*directed){
				// Soergel distance (1 - Jaccard coefficient).
				numerat = neighbFocal.size() + neighbOther.size() - 2*common.size();
				denom = neighbFocal.size() + neighbOther.size() - common.size();
				distm = numerat/denom;

	 			if(*disk){
					row.at(j-i-1) = distm;
				}else{
					dissvec[k] = distm;
					}
			}else if(*weighted && !*directed){
				// Loop through sorted node neighbourhood union and extract corresponding weights.
				for(sit = total.begin(); sit != total.end(); sit++){
					aI.push_back( mapA[*sit] ); // Equals zero if unmatched node key.
					aI.push_back( mapB[*sit] );
					aI.push_back( mapA[*sit] );
					aI.push_back( mapB[*sit] );
					aJ.push_back( mapC[*sit] );
					aJ.push_back( mapD[*sit] );
					aJ.push_back( mapD[*sit] );
					aJ.push_back( mapC[*sit] );
					}

				// Tanimoto coefficient.
				dotprod = inner_product(aI.begin(),aI.end(),aJ.begin(),0.0);
				absA = inner_product(aI.begin(),aI.end(),aI.begin(),0.0);
				absB = inner_product(aJ.begin(),aJ.end(),aJ.begin(),0.0);
				if(*disk){
					row.at(j-i-1) = 1.0 - (dotprod/(absA + absB - dotprod));
				}else{
					dissvec[k] = 1.0 - (dotprod/(absA + absB - dotprod));
					}

			}else if(*weighted && *directed){
				// Calculate directed weights - 1 if in same direction, 0.5 (or dirweight) if not.
				getDirectedWeights_all(dirWeights, common, edgeA, edgeB, edgeA.at(i), edgeA.at(j), numedg, dirweight);
				getDirectedWeights_all(dirWeights, common, edgeA, edgeB, edgeA.at(i), edgeA.at(j), numedg, dirweight);
				
				// Loop through sorted node neighbourhood union and extract corresponding weights.
				for(sit = total.begin(); sit != total.end(); sit++){
					aI.push_back( mapA[*sit] ); // Equals zero if unmatched node key.
					aJ.push_back( mapB[*sit] );
					if(dirWeights[*sit] != 0){
						aI.back() = aI.back()*dirWeights[*sit];
						aJ.back() = aJ.back()*dirWeights[*sit];
						}
					}
				// Tanimoto coefficient.
				dotprod = inner_product(aI.begin(),aI.end(),aJ.begin(),0.0);
				absA = inner_product(aI.begin(),aI.end(),aI.begin(),0.0);
				absB = inner_product(aJ.begin(),aJ.end(),aJ.begin(),0.0);
				if(*disk){
					row.at(j-i-1) = 1.0 - (dotprod/(absA + absB - dotprod));
				}else{
					dissvec[k] = 1.0 - (dotprod/(absA + absB - dotprod));
					}
					
			}else if(!*weighted && *directed){
				// Calculate directed weights - 1 if in same direction, 0.5 if not.
				getDirectedWeights_all(dirWeights, common, edgeA, edgeB, nonshared.at(0), nonshared.at(1), numedg, dirweight);

				// Loop through sorted node neighbourhood union and extract corresponding weights.
				for(sit = total.begin(); sit != total.end(); sit++){
					if(find(neighbFocal.begin(),neighbFocal.end(),*sit) != neighbFocal.end() && dirWeights[*sit] != 0){
						aI.push_back( dirWeights[*sit] );
					}else if(find(neighbFocal.begin(),neighbFocal.end(),*sit) != neighbFocal.end() && dirWeights[*sit] == 0){
						aI.push_back(1);
					}else if(find(neighbFocal.begin(),neighbFocal.end(),*sit) == neighbFocal.end()){
						aI.push_back(0);
						}

					if(find(neighbOther.begin(),neighbOther.end(),*sit) != neighbOther.end() && dirWeights[*sit] != 0){
						aJ.push_back( dirWeights[*sit] );
					}else if(find(neighbOther.begin(),neighbOther.end(),*sit) != neighbOther.end() && dirWeights[*sit] == 0){
						aJ.push_back(1);
					}else if(find(neighbOther.begin(),neighbOther.end(),*sit) == neighbOther.end()){
						aJ.push_back(0);
						}
					}
				// Tanimoto coefficient.
				dotprod = inner_product(aI.begin(),aI.end(),aJ.begin(),0.0);
				absA = inner_product(aI.begin(),aI.end(),aI.begin(),0.0);
				absB = inner_product(aJ.begin(),aJ.end(),aJ.begin(),0.0);
				if(*disk){
					row.at(j-i-1) = 1.0 - (dotprod/(absA + absB - dotprod));
				}else{
					dissvec[k] = 1.0 - (dotprod/(absA + absB - dotprod));
					}
				
				}

			//Rprintf("\ninds.at(j) %d\ncommon.size %d\ntotal.size %d\n",inds.at(j),common.size(),total.size());

			diff.clear();
			nonshared.clear();
			neighbA.clear();
			neighbB.clear();
			neighbOther.clear();
			common.clear();
			total.clear();
			mapA.clear();
			mapB.clear();
			aI.clear();
			aJ.clear();
			dirWeights.clear();

			k = k + 1;

			}

		if(*disk){
			// Collapse strings of '1's in each row to compress file size.
			temp = row;
			sum = 0; runn = 0;
			first = -1; last = -1;

			for(j = 0; j < row.size(); j++){
				if(row.at(j) != 1 && first == -1){
					continue;
				}else if(row.at(j) != 1 && first != -1){
					if(first == last){
						first = -1;
						last = -1;
						sum = 0;
						continue;
						}
					temp.erase(temp.begin() + first - runn, temp.begin() + last + 1 - runn);
					temp.insert(temp.begin() + first - runn, sum);
					runn = runn + sum - 1;
					last = -1;
					first = -1;
					sum = 0;
				}else if(row.at(j) == 1 && last == -1){
					first = j;
					last = j;
					sum++;
				}else if(row.at(j) == 1 && last != -1){
					last++;
					sum++;
				}
				
				if(j == row.size()-1 && last != -1){
					if(first == last){continue;}
					temp.erase(temp.begin() + first - runn, temp.end());
					temp.insert(temp.begin() + first - runn, sum);
					}
				}

			row = temp;

			for(j = 0; j < row.size(); j++){
				outfile << row.at(j) << " ";
				}
			outfile << endl;

			rowlen[i] = row.size();

			}

		row.clear();
		inds.clear();
		neighbFocal.clear();

		}

	if(*disk){
		outfile.close();
		}

	}

	}






