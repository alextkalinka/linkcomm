/* Code to calculate edge similarities for hierarchical clustering of edges to produce link communities.
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


void getNodeNeighbourhood(map<int, set<int> > &nN, map<int, map<int,double> > &wM, vector<int> &eA, vector<int> &eB, int nonSA, int nonSB, unsigned int *numedg, bool *weighted, double *weights, int which)

	{

	unsigned int finda = 0, findab = 0, findb = 0, findba = 0, add;
	bool end = FALSE, first = TRUE;

	if(which == 1){
		findab = *numedg;
		findba = *numedg;
	}else if(which == 2){
		finda = *numedg;
		findb = *numedg;
		}

	// Get first-order node neighbourhoods for the non-shared nodes.
	while(end == FALSE){
		if((finda != *numedg && which == 0) || which == 1){
			if(first){add = 0;}else{add = finda+1;}
			if(add > *numedg){add = *numedg;}
			finda = find(eA.begin()+add, eA.end(), nonSA) - eA.begin();
			if(finda != *numedg){
				nN[nonSA].insert(eB.at(finda));
				if(*weighted){
					wM[nonSA].insert( pair<int,double>(eB.at(finda), weights[finda]) );
					}
				}
			}
		if((findb != *numedg && which == 0) || which == 1){
			if(first){add = 0;}else{add = findb+1;}
			if(add > *numedg){add = *numedg;}
			findb = find(eB.begin()+add, eB.end(), nonSA) - eB.begin();
			if(findb != *numedg){
				nN[nonSA].insert(eA.at(findb));
				if(*weighted){
					wM[nonSA].insert( pair<int,double>(eA.at(findb), weights[findb]) );
					}
				}
			}
		if((findab != *numedg && which == 0) || which == 2){
			if(first){add = 0;}else{add = findab+1;}
			if(add > *numedg){add = *numedg;}
			findab = find(eA.begin()+add, eA.end(), nonSB) - eA.begin();
			if(findab != *numedg){
				nN[nonSB].insert(eB.at(findab));
				if(*weighted){
					wM[nonSB].insert( pair<int,double>(eB.at(findab), weights[findab]) );
					}
				}
			}
		if((findba != *numedg && which == 0) || which == 2){
			if(first){add = 0;}else{add = findba+1;}
			if(add > *numedg){add = *numedg;}
			findba = find(eB.begin()+add, eB.end(), nonSB) - eB.begin();
			if(findba != *numedg){
				nN[nonSB].insert(eA.at(findba));
				if(*weighted){
					wM[nonSB].insert( pair<int,double>(eA.at(findba), weights[findba]) );
					}
				}
			}

		if(finda == *numedg && findb == *numedg && findab == *numedg && findba == *numedg){
			end = TRUE;
		}else{
			first = FALSE;
			}

		}

	}



void getDirectedWeights(map<int,float> &dW, set<int> &comm, vector<int> &eA, vector<int> &eB, int nonSA, int nonSB, unsigned int *numedg, double *dirw)
	
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



void getEdgeSimilarities(int *ea, int *eb, unsigned int *numedg, unsigned int *rowlen, double *weights, bool *directed, double *dirweight, bool *weighted, bool *disk, double *dissvec, bool *bipartite, bool *verbose)

	{

	unsigned int i, j, finda, findab, findb, findba, add, which;
	int first = -1, last = -1, sum = 0, runn = 0;
	double dotprod, absA, absB, numerat, denom, distm;
	float prog;
	bool end = FALSE;
	vector<double> row;
	vector<double> temp;
	vector<int> edgeA;
	vector<int> edgeB;
	vector<int> inds;
	vector<int> nonshared;
	set<int> nodesI;
	set<int> nodesJ;
	set<int> diff;
	set<int> neighbA;
	set<int> neighbB;
	set<int> common;
	set<int> total;
	set<int>::iterator sit;
	map<int, set<int> > nodeNeighb;
	nodeNeighb[-1].insert(0); // Initialize with dummy variable.
	map<int, set<int> >::iterator mitA;
	map<int, set<int> >::iterator mitB;
	map<int, map<int,double> > weightMap;
	weightMap[-1].insert( pair<int,double>(-1, 0.0) ); // Initialize.
	map<int,double> mapA;
	map<int,double> mapB;
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

	// Loop through edges and calculate edge similarities.
	for(i = 0; i < *numedg-1; i++){

		if(*verbose){
			prog = (i+0.0)/(*numedg-2)*100;

			Rprintf("\r   Calculating edge similarities for %d edges... %3.2f%%",*numedg,prog);

			R_FlushConsole();
			R_ProcessEvents();
			}

		row.assign(*numedg-1-i, 1);
		end = FALSE; j = 0;
		finda = 0; findb = 0; findab = 0; findba = 0;
        
		// Find edges that share a node with this edge.
		while(end == FALSE){
			if(finda != *numedg){
				if(j==0){add = i+1;}else{add = finda+1;}
				finda = find(edgeA.begin()+add,edgeA.end(),edgeA.at(i)) - edgeA.begin();
				if(finda != *numedg){
					inds.push_back(finda);
					}
				}
			if(findb != *numedg){
				if(j==0){add = i+1;}else{add = findb+1;}
				findb = find(edgeB.begin()+add,edgeB.end(),edgeB.at(i)) - edgeB.begin();
				if(findb != *numedg){
					inds.push_back(findb);
					}
				}
			if(findab != *numedg){
				if(j==0){add = i+1;}else{add = findab+1;}
				findab = find(edgeB.begin()+add,edgeB.end(),edgeA.at(i)) - edgeB.begin();
				if(findab != *numedg){
					inds.push_back(findab);
					}
				}
			if(findba != *numedg){
				if(j==0){add = i+1;}else{add = findba+1;}
				findba = find(edgeA.begin()+add,edgeA.end(),edgeB.at(i)) - edgeA.begin();
				if(findba != *numedg){
					inds.push_back(findba);
					}
				}
			
			if(finda == *numedg && findb == *numedg && findab == *numedg && findba == *numedg){
				end = TRUE;
			}else{
				j = 1;
				}

			}

		nodesI.insert(edgeA.at(i)); // The nodes for this edge.
		nodesI.insert(edgeB.at(i));

		// Loop through edges that share a node with the current edge and calculate similarity scores.
		for(j = 0; j < inds.size(); j++){

			// Get the two non-shared nodes.
			nodesJ.insert(edgeA.at(inds.at(j)));
			nodesJ.insert(edgeB.at(inds.at(j)));

			//Rprintf("\ninds.size %d\nedgeA.at(inds.at(j)) %d\nedgeB.at(inds.at(j)) %d\n",inds.size(),edgeA.at(inds.at(j)),edgeB.at(inds.at(j)));

			set_difference(nodesI.begin(),nodesI.end(),nodesJ.begin(),nodesJ.end(), inserter(diff, diff.begin()));
			set_difference(nodesJ.begin(),nodesJ.end(),nodesI.begin(),nodesI.end(), inserter(diff, diff.begin()));
			
			for(sit = diff.begin(); sit != diff.end(); sit++){
				nonshared.push_back(*sit);
				}

			//for(k=0;k<nonshared.size();k++){
			//	Rprintf("\nnonsh: %d\n",nonshared.at(k));
			//	}

			mitA = nodeNeighb.find(nonshared.at(0));
			mitB = nodeNeighb.find(nonshared.at(1));

			if(mitA == nodeNeighb.end() && mitB == nodeNeighb.end()){
				which = 0;
				getNodeNeighbourhood(nodeNeighb, weightMap, edgeA, edgeB, nonshared.at(0), nonshared.at(1), numedg, weighted, weights, which);
			}else if(mitA == nodeNeighb.end() && mitB != nodeNeighb.end()){
				which = 1;
				getNodeNeighbourhood(nodeNeighb, weightMap, edgeA, edgeB, nonshared.at(0), nonshared.at(1), numedg, weighted, weights, which);
			}else if(mitA != nodeNeighb.end() && mitB == nodeNeighb.end()){
				which = 2;
				getNodeNeighbourhood(nodeNeighb, weightMap, edgeA, edgeB, nonshared.at(0), nonshared.at(1), numedg, weighted, weights, which);
				}

			neighbA = nodeNeighb[nonshared.at(0)];
			neighbB = nodeNeighb[nonshared.at(1)];

			if(!(*bipartite)){ // Inclusive node neighbourhood only if not bipartite.
				neighbA.insert(nonshared.at(0));
				neighbB.insert(nonshared.at(1));
				}

			if(*weighted){
				mapA = weightMap[nonshared.at(0)];
				mapB = weightMap[nonshared.at(1)];
				if(!(*bipartite)){
					mapA.insert( pair<int,double>(nonshared.at(0), 1) );
					mapB.insert( pair<int,double>(nonshared.at(1), 1) );
					}
				}

			//for(sit = neighbA.begin(); sit != neighbA.end(); sit++){
			//	Rprintf("\nneighbA: %d\n",*sit);
			//	}
			//for(sit = neighbB.begin(); sit != neighbB.end(); sit++){
			//	Rprintf("\nneighbB: %d\n",*sit);
			//	}

			set_intersection(neighbA.begin(),neighbA.end(),neighbB.begin(),neighbB.end(), inserter(common, common.begin()));

			if(*weighted || *directed){
				set_union(neighbA.begin(),neighbA.end(),neighbB.begin(),neighbB.end(), inserter(total, total.begin()));
				}

			if(!*weighted && !*directed){
				// Soergel distance (1 - Jaccard coefficient).
				numerat = neighbA.size() + neighbB.size() - 2*common.size();
				denom = neighbA.size() + neighbB.size() - common.size();
				distm = numerat/denom;

				if(*disk){
					row.at(inds.at(j)-i-1) = distm;
				}else{
					dissvec[inds.at(j)-i-1+runn] = distm;
					}
			}else if(*weighted && !*directed){
				// Loop through sorted node neighbourhood union and extract corresponding weights.
				for(sit = total.begin(); sit != total.end(); sit++){
					aI.push_back( mapA[*sit] ); // Equals zero if unmatched node key.
					aJ.push_back( mapB[*sit] );
					}

				// Tanimoto coefficient.
				dotprod = inner_product(aI.begin(),aI.end(),aJ.begin(),0.0);
				absA = inner_product(aI.begin(),aI.end(),aI.begin(),0.0);
				absB = inner_product(aJ.begin(),aJ.end(),aJ.begin(),0.0);
				if(*disk){
					row.at(inds.at(j)-i-1) = 1.0 - (dotprod/(absA + absB - dotprod));
				}else{
					dissvec[inds.at(j)-i-1+runn] = 1.0 - (dotprod/(absA + absB - dotprod));
					}

			}else if(*weighted && *directed){
				// Calculate directed weights - 1 if in same direction, 0.5 if not.
				getDirectedWeights(dirWeights, common, edgeA, edgeB, nonshared.at(0), nonshared.at(1), numedg, dirweight);
				
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
					row.at(inds.at(j)-i-1) = 1.0 - (dotprod/(absA + absB - dotprod));
				}else{
					dissvec[inds.at(j)-i-1+runn] = 1.0 - (dotprod/(absA + absB - dotprod));
					}
					
			}else if(!*weighted && *directed){
				// Calculate directed weights - 1 if in same direction, 0.5 if not.
				getDirectedWeights(dirWeights, common, edgeA, edgeB, nonshared.at(0), nonshared.at(1), numedg, dirweight);

				// Loop through sorted node neighbourhood union and extract corresponding weights.
				for(sit = total.begin(); sit != total.end(); sit++){
					if(find(neighbA.begin(),neighbA.end(),*sit) != neighbA.end() && dirWeights[*sit] != 0){
						aI.push_back( dirWeights[*sit] );
					}else if(find(neighbA.begin(),neighbA.end(),*sit) != neighbA.end() && dirWeights[*sit] == 0){
						aI.push_back(1);
					}else if(find(neighbA.begin(),neighbA.end(),*sit) == neighbA.end()){
						aI.push_back(0);
						}

					if(find(neighbB.begin(),neighbB.end(),*sit) != neighbB.end() && dirWeights[*sit] != 0){
						aJ.push_back( dirWeights[*sit] );
					}else if(find(neighbB.begin(),neighbB.end(),*sit) != neighbB.end() && dirWeights[*sit] == 0){
						aJ.push_back(1);
					}else if(find(neighbB.begin(),neighbB.end(),*sit) == neighbB.end()){
						aJ.push_back(0);
						}
					}
				// Tanimoto coefficient.
				dotprod = inner_product(aI.begin(),aI.end(),aJ.begin(),0.0);
				absA = inner_product(aI.begin(),aI.end(),aI.begin(),0.0);
				absB = inner_product(aJ.begin(),aJ.end(),aJ.begin(),0.0);
				if(*disk){
					row.at(inds.at(j)-i-1) = 1.0 - (dotprod/(absA + absB - dotprod));
				}else{
					dissvec[inds.at(j)-i-1+runn] = 1.0 - (dotprod/(absA + absB - dotprod));
					}
				
				}

			//Rprintf("\ninds.at(j) %d\ncommon.size %d\ntotal.size %d\n",inds.at(j),common.size(),total.size());

			nodesJ.clear();
			diff.clear();
			nonshared.clear();
			neighbA.clear();
			neighbB.clear();
			common.clear();
			total.clear();
			mapA.clear();
			mapB.clear();
			aI.clear();
			aJ.clear();
			dirWeights.clear();

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

		}else{
			runn = runn + *numedg - 1 - i;			
			}

		row.clear();
		inds.clear();
		nodesI.clear();

		}

	if(*disk){
		outfile.close();
		}

	}

	}






