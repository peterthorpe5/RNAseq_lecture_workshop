#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>
#include <map>
#include <set>
#include <cstdlib>

using namespace std;
ostringstream outfile, outfilep;

int main (int argc, char **argv)
{

	if(argc != 6)
	{
		cerr << "\nSyntax: select_results RNentropyfile.summary.res GPV_threshold LPV_threshold sample_num rep_num" << endl << endl;
		exit(EXIT_FAILURE);
	}	

	multimap<double, string> LINES, LINES_BENJAMINI;
	vector<string> labels;
	unsigned int sam_num, rep_num;
	double G_T, L_T;
	string line, t;

	G_T = atof(argv[2]);
        G_T = -log10(G_T);
	L_T = atof(argv[3]);
        L_T = -log10(L_T);
	sam_num = atoi(argv[4]);
	rep_num = atoi(argv[5]);

	for(unsigned int i = 0; i < sam_num; i++)
	{
		ostringstream sstr;

		sstr << "SAMPLE_" << i+1;
		labels.push_back(sstr.str());
	}

	ifstream in(argv[1]);

	if(!in)
	{
		cerr << "\nCan't find " << argv[1] << endl;
		exit(EXIT_FAILURE);
	}

	outfile << argv[1] << ".selected";	
	outfilep << argv[1] << ".pmi";

	getline(in,line);

	istringstream buf(line);
	unsigned int lcount = 0, ccount = 0;
	bool gflag = false;

	while(buf >> t)
	{
		if(t != "GL_LPV")
			lcount++;
		else
		{
			gflag = true;
			break;
		}

		if(t.find("COMMENT_") != string::npos)
			ccount++;
	}

	if(lcount != 3 + ccount)
	{
		cerr << endl << endl << argv[1] << " does not look like a valid RNentropy summary file\n" << endl;
		exit(EXIT_FAILURE);
	}

	while(getline(in,line))
	{
		if(line.empty())
			continue;

		if(line[0] == '#')
			continue;

		istringstream str(line);
		string trash;
		double GPV;

		for(unsigned int i = 0; i < lcount; i++)
			str >> trash;

		str >> GPV;

		LINES.insert(make_pair(pow(10, -GPV), line));
		
	}	

	in.close();
	
//	cerr << "\nLINES.size() = " << LINES.size() << endl;

	multimap<double, string>::iterator mi = LINES.begin();
	unsigned int count = 1;
	ofstream sout(outfile.str().c_str(), ios::out), pout(outfilep.str().c_str(), ios::out);

	while(mi != LINES.end())
	{
		double benj_G_PV = mi->first * ((double)LINES.size() / (double)count);

		if(benj_G_PV > 1)
			benj_G_PV = 1;
		
		LINES_BENJAMINI.insert(make_pair(benj_G_PV, mi->second));

//		sout << mi->first << '\t' << benj_G_PV << '\t' << mi->second << endl;
	
		mi++;
		count++;
	}		

	mi =  LINES_BENJAMINI.begin();	

	sout << "GENE_ID\tTR_ID\t---\tGL_LPV\tCORR_GL_LPV\t---\t" ; 

	for(unsigned int i = 0; i < labels.size(); i++)
		sout << labels.at(i) << '\t';

	sout << endl;

	vector<vector<unsigned int> > bv(labels.size(), vector<unsigned int>());

	while(mi != LINES_BENJAMINI.end())
	{
		istringstream str(mi->second);

		double Log_benj_G_PV = -log10(mi->first);		

		if(Log_benj_G_PV < G_T)
			break;

		string gene_id, tr_id, trash;
		double or_pv;

		str >> gene_id >> tr_id >> trash;

		for(unsigned int i = 0; i < ccount; i++)	
			str >> trash;

		str >> or_pv >> trash;

		sout << gene_id << '\t' << tr_id << '\t' << "---" << '\t' << or_pv << '\t' << Log_benj_G_PV << '\t' << "---" << '\t';
		unsigned int scount = 0;

		while(str)
		{
			vector<int> flag;			

			for(unsigned int i = 0; i < rep_num; i++)
			{	
				double lpv;
				str >> lpv;

				if(fabs(lpv) < L_T)
					flag.push_back(0);
				else
					if(lpv > 0)
						flag.push_back(1);
					else
						flag.push_back(-1);
						
			}	

			if(str)
			{
				int fsum = 0, znum = 0;

				for(unsigned int i = 0; i < rep_num; i++)
				{
					fsum += flag[i];
				
					if(flag[i] == 0)
						znum++;
				}

				if(fsum == rep_num)
				{
					sout << 1 << '\t';
					bv[scount].push_back(1);
				}
				else if(fsum == -(int)rep_num)
				{
					sout << -1 << '\t';
					bv[scount].push_back(0);
				}
				else if(znum > 0)
				{
					sout << 0 << '\t';
					bv[scount].push_back(0);
				}
				else
				{
					sout << 'X' << '\t';
					bv[scount].push_back(0);
				}
			}

			scount++;
		}	

		sout << endl;
	
		mi++;
	}

	sout.close();
	
	double **pmi = new double*[bv.size()];	
	double **npmi = new double*[bv.size()];

	for(unsigned int i = 0; i < bv.size(); i++)
	{
		pmi[i] = new double[bv.size()];
		npmi[i] = new double[bv.size()];
	}

	for(unsigned int x = 0; x < bv.size(); x++)
	{
		double xsum = 0, xfreq;

		for(unsigned int i = 0; i < bv[x].size(); i++)
			xsum += bv[x][i];

		xfreq = xsum / (double)bv[x].size();

		for(unsigned int y = 0; y < bv.size(); y++)
		{
			if(x > y)
			{
				pmi[x][y] = pmi[y][x];
				npmi[x][y] = npmi[y][x];
				continue;
			}

			double ysum = 0, yfreq, xysum = 0, xyfreq, h_xy;

                	for(unsigned int i = 0; i < bv[y].size(); i++)
			{
                        	ysum += bv[y][i];
				
				xysum += (bv[x][i] * bv[y][i]);
			}

                	yfreq = ysum / (double)bv[y].size();
			xyfreq = xysum / (double)bv[y].size();

			h_xy = log2(1/xyfreq);

			pmi[x][y] = log2(xyfreq / (xfreq * yfreq));
			npmi[x][y] = pmi[x][y] / h_xy;
		}
	}

	pout << "***POINT MUTUAL INFORMATION TABLE" << endl;

	pout << '\t';

	for(unsigned int i = 0; i < labels.size(); i++)
		pout << labels.at(i).substr(8) << '\t';

	pout << endl;

	for(unsigned int y = 0; y < bv.size(); y++)
	{
		pout << labels.at(y).substr(8) << '\t';

		for(unsigned int x = 0; x < bv.size(); x++)
			pout << pmi[x][y] << '\t';

		pout << endl;
	}

	pout << endl << endl << endl;

	pout << "***NORM. POINT MUTUAL INFORMATION TABLE" << endl;

        pout << '\t';

        for(unsigned int i = 0; i < labels.size(); i++)
                pout << labels.at(i).substr(8) << '\t';

        pout << endl;

        for(unsigned int y = 0; y < bv.size(); y++)
        {
                pout << labels.at(y).substr(8) << '\t';

                for(unsigned int x = 0; x < bv.size(); x++)
                        pout << npmi[x][y] << '\t';

                pout << endl;
        }

	
	for(unsigned int i = 0; i < bv.size(); i++)
	{
                delete[] pmi[i];
		delete[] npmi[i];
	}

	delete [] pmi;
	delete [] npmi;

	return EXIT_SUCCESS;
}
