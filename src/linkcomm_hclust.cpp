/* Code to carry out single linkage hierarchical clustering of edges in a network.
 * Will be loaded as a shared object into R (dyn.load('name.so')) and called from within R using .C("hclustLinkComm",...).
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
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <sstream>
#include <map>

#include <R.h>

extern "C" {

using namespace std;


void compressRow(vector<float> &v)

	{
	// Compresses 1's in a row of the diss matrix.

	unsigned int i;
	int runn = 0, add = 0, first = -1, last = -1;
	float sum = 0;
	vector<float> temp = v;

	for(i = 0; i < v.size(); i++){
		if(v.at(i) < 1 && first == -1){
			continue;
		}else if(v.at(i) < 1 && first != -1){
			if(first == last){
				first = -1;
				last = -1;
				sum = 0;
				add = 0;
				continue;
				}
			temp.erase(temp.begin() + first - runn, temp.begin() + last + 1 - runn);
			temp.insert(temp.begin() + first - runn, sum);
			runn = runn + add - 1;
			last = -1;
			first = -1;
			sum = 0;
			add = 0;
		}else if(v.at(i) >= 1 && last == -1){
			first = i;
			last = i;
			sum = sum + v.at(i);
			add++;
		}else if(v.at(i) >= 1 && last != -1){
			last++;
			sum = sum + v.at(i);
			add++;
			}
			
		if(i == v.size()-1 && last != -1){
			if(first == last){continue;}
			temp.erase(temp.begin() + first - runn, temp.end());
			temp.insert(temp.begin() + first - runn, sum);
			}
		}
	
	v = temp;
	temp.clear();

	}


void hclustLinkComm(unsigned int *numedg, int *rowlen, float *heights, int *hca, int *hcb, bool *verbose)
	{
	
	int row, col, numedgU = *numedg;
	float diff, best, min, prog, *arr;
	vector<int> merges; // Edge pairs that are to be agglomerated.
	vector<int> einds; // Edge indices to be updated as edges are deleted and new clusters added.
	vector<int>::iterator eit; // Vector iterators.
	vector<int>::iterator frit;
	vector<int>::iterator fcit; 
	vector<int> rowL; // Length of rows to be updated as they change and are deleted.
	vector<int> tempL;
	map<int, map<int,float> > oldrow;
	map<int, map<int,float> > oldcol;
	map<int,float>::iterator mrit;
	map<int,float>::iterator mcit;
	vector<float> rowM;
	vector<float> endM;
	int wleft, fin = 0;
	unsigned int i, j, m, count = 0, numM = 0;
	int k, p;
	bool sones = FALSE;
	string line;

	ifstream infile;
	ofstream outfile;
	stringstream ss (stringstream::in | stringstream::out);

	for(i = 0; i < *numedg; i++){ // Set up edge indices.
		einds.push_back( -1*(i+1) );
		}

	copy(rowlen, rowlen + (*numedg-1), back_inserter(rowL)); // Copy row lengths into a vector.

	infile.open("linkcomm_diss.txt", ios::in );
	if(! infile.is_open()){
			Rprintf("\nERROR: linkcomm_diss.txt not found!\n"); return;
			}

	//remove("linkcomm_diss");
	outfile.open("linkcomm_diss", ios::out | ios::binary);

	/* Read in decimal data file and write out binary data file. */
	while(getline(infile,line)){
		ss << line;
		while(ss >> diff){
			outfile.write((char *) &diff, sizeof(float));
			}
		ss.clear();
		}

	infile.close();	
	remove("linkcomm_diss.txt");

	outfile.close();

	/* There are (numedg-1) agglomerations carried out. 
	 * Ties in independent pairs are merged simultaneously. 
	 * Ties in the same pairs are agglomerated according to their order in the diss matrix. */
	while(numM < *numedg-1){

		if(*verbose){		
			prog = (numM+0.0)/(*numedg-2)*100;

			Rprintf("\r   Hierarchical clustering of edges... %3.2f%%",prog);

			R_FlushConsole();
			R_ProcessEvents();
			}

		best = 2.0; k = 0;
		row = 0; col = 0;

		infile.open("linkcomm_diss", ios::in | ios::binary);
		if(! infile.is_open()){
			Rprintf("\nERROR: linkcomm_diss not found!\n"); return;
			}

		/* Two loops - one reads data and finds clusters, the other writes the updated diss matrix to file. */
		while(row < numedgU-1){
			arr = new float[rowL.at(row)];
			infile.read((char *) arr, rowL.at(row)*sizeof(float)); // Read in row.
			while(col <= numedgU-2){ // Handle row.
				col++;
				diff = arr[k]; 
				if(diff > 1){ // Handle summed '1's.
					diff = 1; sones = TRUE;
					}
				if(diff < best){ // Update best.
					best = diff;
					merges.erase(merges.begin(),merges.end());
					merges.push_back(row);
					merges.push_back(col);
				}else if(diff == best){
					frit = find(merges.begin(),merges.end(),row); // Is this row and column already merged?
					fcit = find(merges.begin(),merges.end(),col);
					if(frit == merges.end() && fcit == merges.end()){
						merges.push_back(row);
						merges.push_back(col);
						}
					}
				if(sones == TRUE){
					col = col + arr[k] - 1;
					sones = FALSE;
					}
				k++;
				}
			row++;
			col = row; k = 0;
			delete [] arr;
			arr = NULL;
			}
		//Rprintf("\ncol is %d, row is %d, and numM is %d \n",col,row,numM);
		//Rprintf("there are %d merges\n",merges.size());
		//for(i=0;i<merges.size();i++){
		//	Rprintf("%d\n",merges.at(i));
		//	}
		//if(numM==38){return;}
		//Rprintf("\n%1.3f",best);

		infile.clear();
		infile.seekg(0, ios::beg);
		outfile.open("temp", ios::out | ios::binary);
		if(! outfile.is_open()){
			Rprintf("\nERROR: cannot open temp!\n"); return;
			}

		row = 0; col = 0;
		tempL.assign(rowL.size(),0);

		wleft = numedgU - merges.size() + (merges.size()/2) - 1; // Number of new entries to be written.

		for(i = 0; i < (merges.size()/2); i++){ // Set up empty old row and col maps of maps.
			oldrow.insert( pair<int,map<int,float> >() );
			oldcol.insert( pair<int,map<int,float> >() );
			}


		while(row < numedgU-1){
			arr = new float[rowL.at(row)];
			infile.read((char *) arr, rowL.at(row)*sizeof(float)); // Read in row.
			while(col <= numedgU-2){ // Handle row.
				col++;
				diff = arr[k];
				frit = find(merges.begin(),merges.end(),row); // Is this row and column merged?
				fcit = find(merges.begin(),merges.end(),col);
				if(diff > 1){ // Handle summed '1's.
					//p = 0;
					for(i = 0; i < merges.size(); i++){
						//if(i % 2 != 0){p++;}
						if(col <= merges.at(i) && merges.at(i) < (col + arr[k])){
							if(i % 2 == 0){
								//oldrow[(i/2)].insert( pair<int,float>(row,1) );
								count++;
							}else{
								if(row != merges.at(i-1)){
									//oldcol[(i-p)].insert( pair<int,float>(row,1) );
									count++;
									}
								}
							}
						}
					diff = diff - count; // Update length of summed '1's.
					if(diff == 0){ // Row length shortens by 1.
						//tempL.at(row)--;
						diff = -1;
					}else if(diff == 1 && frit == merges.end()){
						//outfile.write((char *) &diff, sizeof(float));
						rowM.push_back(diff);
						diff = -1;
					}else if(frit == merges.end()){
						//outfile.write((char *) &diff, sizeof(float));
						rowM.push_back(diff);
						}
					col = col + arr[k] - 1;
					count = 0;					
					}
				if(frit == merges.end() && fcit == merges.end() && diff <= 1 && diff != -1){ // Save non-clustered diff.
					//outfile.write((char *) &diff, sizeof(float));
					rowM.push_back(diff);
					//Rprintf("diff: %1.3f\n",diff);
				}else if((frit != merges.end() && diff <= 1 && diff != -1) || (fcit != merges.end() && diff <= 1 && diff != -1)){
					p = 0;
					for(i = 0; i < merges.size(); i++){ // Save clustered diff.
						if(i % 2 != 0){p++;}
						if(row == merges.at(i)){
							if(i % 2 == 0){
								oldrow[(i/2)].insert( pair<int,float>(col,diff) );
								//tempL.at(row)--;
								//Rprintf("diff: %1.3f\n",diff);
							}else{
								oldcol[(i-p)].insert( pair<int,float>(col,diff) );
								//tempL.at(row)--;
								//Rprintf("diff: %1.3f\n",diff);
								}
						}else if(col == merges.at(i)){
							if(i % 2 == 0){
								oldrow[(i/2)].insert( pair<int,float>(row,diff) );
								//tempL.at(row)--;
								//Rprintf("diff: %1.3f\n",diff);
							}else{
								oldcol[(i-p)].insert( pair<int,float>(row,diff) );
								//tempL.at(row)--;
								//Rprintf("diff: %1.3f\n",diff);
								}
							}
						}
					}
				if(col == numedgU-1 && row != numedgU-2){ // End of row; write row with new cluster minima.
					if(frit == merges.end()){ // Delete merged rows.

						//Rprintf("row: %d, col %d, oldind %d, oldrow.size %d, oldcol.size %d\n",row,col,oldind,oldrow.size(),oldcol.size());
						for(i = 0; i < (merges.size()/2); i++){

							mrit = oldrow[i].find(row);
							mcit = oldcol[i].find(row);
							if(mrit == oldrow[i].end() && mcit == oldcol[i].end()){ // Not there, so both are 1.
								rowM.push_back(1);
							}else if(mrit != oldrow[i].end() && mcit == oldcol[i].end()){
								if(oldrow[i][row] >= 1){
									rowM.push_back(1);
								}else{
									rowM.push_back(oldrow[i][row]);
									}
							}else if(mrit == oldrow[i].end() && mcit != oldcol[i].end()){
								if(oldcol[i][row] >= 1){
									rowM.push_back(1);
								}else{
									rowM.push_back(oldcol[i][row]);
									}
							}else if(oldrow[i][row] < oldcol[i][row]){
								if(oldrow[i][row] >= 1){
									rowM.push_back(1);
								}else{
									rowM.push_back(oldrow[i][row]);
									}
							}else{
								if(oldcol[i][row] >= 1){
									rowM.push_back(1);
								}else{
									rowM.push_back(oldcol[i][row]);
									}
								}

							if(mrit != oldrow[i].end()){ // Flush out data.
								oldrow[i].erase(row);
								}
							if(mcit != oldcol[i].end()){
								oldcol[i].erase(row);
								}

							}

						compressRow(rowM);

						// Write out row in one go.
						for(i = 0; i < rowM.size(); i++){
							min = rowM.at(i);
							outfile.write((char *) &min, sizeof(float));
							tempL.at(row)++;
							//Rprintf("min %1.3f\n",min);
							}

						rowM.clear();

						wleft--;
						}
					}
				if(row == numedgU-2){ // Handle last diff in file.
					count = wleft - (merges.size()/2) + 1;
					//Rprintf("wleft: %d, count = %d\n",wleft,count);
					for(i = count; i > 0; i--){ // Write remaining "non-merged" minima.
						for(j = 0; j < (merges.size()/2); j++){
							if(find(merges.begin(),merges.end(),numedgU-1) != merges.end()){ // Correct index.
								p = numedgU - i - 1;
							}else{
								p = numedgU - i;
								}
							
							mrit = oldrow[j].find(p);
							mcit = oldcol[j].find(p);
							if(mrit == oldrow[j].end() && mcit == oldcol[j].end()){ // Not there, so both are 1.
								rowM.push_back(1);
							}else if(mrit != oldrow[j].end() && mcit == oldcol[j].end()){
								if(oldrow[j][p] >= 1){
									rowM.push_back(1);
								}else{
									rowM.push_back(oldrow[j][p]);
									}
							}else if(mrit == oldrow[j].end() && mcit != oldcol[j].end()){
								if(oldcol[j][p] >= 1){
									rowM.push_back(1);
								}else{
									rowM.push_back(oldcol[j][p]);
									}
							}else if(oldrow[j][p] < oldcol[j][p]){
								if(oldrow[j][p] >= 1){
									rowM.push_back(1);
								}else{
									rowM.push_back(oldrow[j][p]);
									}
							}else{
								if(oldcol[j][p] >= 1){
									rowM.push_back(1);
								}else{
									rowM.push_back(oldcol[j][p]);
									}
								}

							if(mrit != oldrow[j].end()){ // Flush out data.
								oldrow[j].erase(p);
								}
							if(mcit != oldcol[j].end()){
								oldcol[j].erase(p);
								}

							}

						compressRow(rowM);

						// Write.
						for(m = 0; m < rowM.size(); m++){
							min = rowM.at(m);
							outfile.write((char *) &min, sizeof(float));
							fin++;
							}

						tempL.push_back(fin);
						fin = 0;
						rowM.clear();		
						p = 0;

						}

					count = (merges.size()/2) - 1;
					for(i = 1; i <= count; i++){ // Write "merged" minima.						
						for(j = i; j <= count; j++){
							for(m = 2*j; m < (2*j+2); m++){ // Assemble merged pairs into vector.
								//Rprintf("%d\n%d\n%d\n",i,j,m);

								mrit = oldrow[i-1].find(merges.at(m));
								mcit = oldcol[i-1].find(merges.at(m));
								if(mrit == oldrow[i-1].end() && mcit == oldcol[i-1].end()){ // Not there, so both are 1.
									rowM.push_back(1);
									rowM.push_back(1);
								}else if(mrit != oldrow[i-1].end() && mcit == oldcol[i-1].end()){
									rowM.push_back(oldrow[i-1][merges.at(m)]);
									rowM.push_back(1);
								}else if(mrit == oldrow[i-1].end() && mcit != oldcol[i-1].end()){
									rowM.push_back(oldcol[i-1][merges.at(m)]);
									rowM.push_back(1);
								}else{
									rowM.push_back(oldrow[i-1][merges.at(m)]);
									rowM.push_back(oldcol[i-1][merges.at(m)]);
									}

								}

							//min = *min_element(endM.begin(),endM.end());

							min = *min_element(rowM.begin(),rowM.end());
							if(min > 1){ min = 1;}
							endM.push_back(min);

							//Rprintf("min %1.3f\n",min);

							//outfile.write((char *) &min, sizeof(float));
							//Rprintf("min: %1.3f\n",min);
							rowM.clear();

							}

						compressRow(endM);

						// Write.
						for(m = 0; m < endM.size(); m++){
							min = endM.at(m);
							outfile.write((char *) &min, sizeof(float));
							fin++;
							}

						tempL.push_back(fin);
						fin = 0;
						endM.clear();
						
						}

					count = 0;

					}
					
				k++;
				//Rprintf("%d ",col);
				}
			row++;
			rowM.clear();
			//if(row==21){return;}
			col = row; k = 0;
			delete [] arr;
			arr = NULL;
			}

		//Rprintf("\ncol is %d, row is %d, and numM is %d \n",col,row,numM);
		//Rprintf("\n");
		//for(i=0;i<oldcol.size();i++){
		//	Rprintf("%1.3f,",oldcol.at(i));
		//	}
		//if(numM==18){
		//	for(i=0;i<oldrow[0].size();i++){
		//		Rprintf("oldrow: %1.3f, ",oldrow[0].at(i));
		//		Rprintf("oldcol: %1.3f\n",oldcol[0].at(i));
		//		}
		//	}

		rowL = tempL;
		tempL.clear();

		//Rprintf("\n");
		//for(i=0;i<rowL.size();i++){
		//	Rprintf("%d,",rowL.at(i));
		//	}
		//Rprintf("\n");

		i = 0; j = 0;
		while(i < merges.size()){
			hca[numM+j] = einds.at(merges.at(i)); // Cluster IDs for this loop (rows always in left-hand column).
			hcb[numM+j] = einds.at(merges.at(i+1));
			i = i + 2;
			j++;
			}

		sort(merges.begin(),merges.end()); // Sort merges for deletion.

		for(i = 0; i < merges.size()/2; i++){ // Add new cluster(s) and heights.
			einds.push_back(numM+i+1);
			heights[numM+i] = best;
			}

		for(i = 0; i < merges.size(); i++){ // Delete merged edge indices and row lengths.
			eit = einds.begin() + merges.at(i) - i;
			//Rprintf("erase: %d\n",einds.at(merges.at(i)-i));
			einds.erase(eit);
			if(merges.at(i) != numedgU-1){
				eit = rowL.begin() + merges.at(i) - i;
				rowL.erase(eit);
				}
			}

		tempL = rowL;
		j = 0;
		for(i = 0; i < rowL.size(); i++){ // Delete zero-length rows.
			if(rowL.at(i) == 0){
				eit = tempL.begin() + i - j;
				tempL.erase(eit);
				j++;
				}
			}

		rowL = tempL;
		tempL.clear();

		numM = numM + (merges.size()/2); // Updated number of merges.
		numedgU = *numedg - 2*numM + numM; // Updated number of edges after agglomeration.
		
		//Rprintf("\n");
		//for(i=0;i<rowL.size();i++){
		//	Rprintf("%d,",rowL.at(i));
		//	}
		//Rprintf("size: %d\n",rowL.size());
		
		merges.clear(); // Clear out old data.
		oldrow.clear();
		oldcol.clear();

		infile.close();
		remove("linkcomm_diss");

		outfile.close();

		if(rename("temp", "linkcomm_diss") != 0){
			Rprintf("\nERROR: temp could not be renamed to linkcomm_diss!\n"); return;
			}


		}

	remove("linkcomm_diss");

	}

	}


