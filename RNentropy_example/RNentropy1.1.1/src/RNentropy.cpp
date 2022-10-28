#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <math.h>
#include <map>
#include <set>
#include <cstdlib>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

using namespace std;

enum 
{
	C_UNKNOWN,
	C_GENE_COL,
	C_TR_COL,
	C_COMMENTS,
	C_EXP,
	C_END	
};

const int BIGLOG = 300;
const double PSEUDOCOUNT = 0.01;
int TR_COL = 0, GENE_COL = 0;
set<unsigned int> S_COL;

struct data_row
{
	string lrow, tr_id, gene_id;
	bool valid;
	vector<string> v_label;
	vector<double> V;
};

void get_data_from_file(string *, vector<bool> *, vector<bool> *, vector<data_row> *, map<string, vector<unsigned int> > *, map<unsigned int, unsigned int>*);
//void write_output_iso(vector<data_row> *, map<string, vector<unsigned int> > *, map<unsigned int, string> *, map<string, set<unsigned int> > *, map<unsigned int, unsigned int> *); //ISO FREQ GTEST
void write_output(vector<data_row> *, map<unsigned int, string> *, map<string, set<unsigned int> > *, map<unsigned int, unsigned int> *); //GLOBAL AND LOCAL SPECIFICITY TEST;
void parse_row(string *, vector<bool>*, vector<bool>*, data_row *, map<unsigned int, unsigned int> *);
void command_line_parser(int, char**);
//void calculate_pv(vector<data_row> *);
void col_num_parser(vector<bool>*);
void col_lab_parser(vector<bool>*);
void get_header_from_file(map<string, set<unsigned int> >*, map<unsigned int, string> *);
unsigned int command_interpreter(string);

string COL, LCOL = "", DATA_FILE;

int main(int argc, char **argv)
{
	vector<data_row> DATA;
	map<string, set<unsigned int> > EXP_REP;
	map<unsigned int, string> COL_TO_EXP;
	map<unsigned int, unsigned int> ARRPOS_TO_COL;
	map<string, vector<unsigned int> > GENE_MAP;

	vector<bool> VB_COL, VB_COL_LAB;
        command_line_parser(argc, argv);

	get_header_from_file(&EXP_REP, &COL_TO_EXP);

	if(!TR_COL)
	{
		cerr << endl << "ERROR: Transcript names column is missing...\n" << endl;
		exit(EXIT_FAILURE);
	}	
	if(!GENE_COL)
        {         
                cerr << endl << "ERROR: Gene names column is missing...\n" << endl;
                exit(EXIT_FAILURE);
        }

        col_num_parser(&VB_COL);
        col_lab_parser(&VB_COL_LAB);	

	if(VB_COL.size() < VB_COL_LAB.size())
                VB_COL.resize(VB_COL_LAB.size(), false);
        else if(VB_COL.size() > VB_COL_LAB.size())
                VB_COL_LAB.resize(VB_COL.size(), false);

	get_data_from_file(&DATA_FILE, &VB_COL, &VB_COL_LAB, &DATA, &GENE_MAP, &ARRPOS_TO_COL);
	
//	calculate_pv(&DATA);

//	write_output(&DATA, &GENE_MAP, &COL_TO_EXP, &EXP_REP, &ARRPOS_TO_COL);	
	write_output(&DATA, &COL_TO_EXP, &EXP_REP, &ARRPOS_TO_COL); 

	return EXIT_SUCCESS;
}

void get_header_from_file(map<string, set<unsigned int> > *EXP_REP, map<unsigned int, string> *COL_TO_EXP)
{
	ifstream in(DATA_FILE.c_str());
	string line;

	while(getline(in,line))
	{
		if(line.empty())
			continue;

		if(line[0] != '#')
			continue;
		else
		{
			istringstream str(line);

			string trash, command, buf1, buf2;

			str >> trash >> command;

			switch(command_interpreter(command))
			{
				case C_END:
					return;
					break;
				case C_GENE_COL:
					str >> buf1;	
					GENE_COL = atoi(buf1.c_str());
					break;
				case C_TR_COL:
					str >> buf1;
                                        TR_COL = atoi(buf1.c_str());
                                        break;
				case C_COMMENTS:
					str >> buf1;
					LCOL = buf1;	
					break;
				case C_EXP:
				{
					str >> buf1 >> buf2;
					map<string, set<unsigned int> >::iterator mi = EXP_REP->find(buf1);

					if(mi == EXP_REP->end())
						mi = EXP_REP->insert(make_pair(buf1, set<unsigned int>())).first;

					istringstream str2(buf2);				
					unsigned int buf;

					while(str2 >> buf)
					{
						if(buf == 0)
						{
							cerr << "\n\nError in header line: " << endl << line << endl << endl;
							cerr << "Column numbering is 1 based" << endl;
							exit(EXIT_FAILURE);
						}

						mi->second.insert(buf);

						S_COL.insert(buf);
		
						COL_TO_EXP->insert(make_pair(buf, buf1));

						if(str2.peek() == ',')
							str2.ignore();	
					}
	
					break;		
				}

				case C_UNKNOWN:
					cerr << endl << "\nUnrecognized line:\n" << line << endl << endl;
					exit(EXIT_FAILURE);
					break;
			}
		}

	}

	in.close();

	return;
}

unsigned int command_interpreter(string command)
{
	if(command == "END")
		return C_END;

	if(command == "GENE_COL")
		return C_GENE_COL;

	if(command == "TR_COL")
		return C_TR_COL;
	
	if(command == "COMMENTS")
		return C_COMMENTS;

	if(command == "EXP")
		return C_EXP;
	

	return C_UNKNOWN;
};

/*void write_output(vector<data_row> *DATA)
{
	for(unsigned int i = 0; i < DATA->size(); i++)
	{
		for(unsigned int j = 0; j < DATA->at(i).v_label.size(); j++)
			cout << DATA->at(i).v_label.at(j) << '\t';

		for(unsigned int j = 0; j < DATA->at(i).V.size(); j++)
                        cout << DATA->at(i).V.at(j) << '\t';

		cout << "\t---\t";

		if(DATA->at(i).mu != 0)
		{
			if(DATA->at(i).pv != 0)
				cout << -log10(DATA->at(i).pv);
			else
				cout << BIG_LOG;
		}
		else 
			cout << 0;

		cout << "\t---\t";

		if(DATA->at(i).mu != 0)
			for(unsigned int j = 0; j < DATA->at(i).VP.size(); j++)
			{
				if(DATA->at(i).V.at(j) > DATA->at(i).mu)
					cout << -log10(DATA->at(i).VP.at(j)) << '\t';
				else
					cout << log10(DATA->at(i).VP.at(j)) << '\t';
			}
		else
			for(unsigned int j = 0; j < DATA->at(i).VP.size(); j++)
				cout << 0 << '\t';

		cout << "\t---\t" << DATA->at(i).mu << endl;
	}

	return;
}
*/
void write_output(vector<data_row> *DATA, map<unsigned int, string> *COL_TO_EXP, map<string, set<unsigned int> > *EXP_REP, map<unsigned int, unsigned int> *ARRPOS_TO_COL)
{
	ofstream out1, out2;
        string fout1, fout2;

        fout1 = fout2 = DATA_FILE;

        fout1 += ".main.res";
        fout2 += ".summary.res";

        out1.open(fout1.c_str());
        out2.open(fout2.c_str());
	
	out1 << "GENE_ID" << '\t' << "TR_ID" << '\t';
	out2 << "GENE_ID" << '\t' << "TR_ID" << '\t';

	for(unsigned int j = 0; j < DATA->at(0).v_label.size(); j++)
	{
		out1 << "COMMENT_" << j+1 << '\t';	
		out2 << "COMMENT_" << j+1 << '\t';
	}

	string pr_exp = "";
	unsigned int rep_c = 0;

	for(unsigned int j = 0; j < DATA->at(0).V.size(); j++)	
	{
		if(pr_exp != COL_TO_EXP->at(ARRPOS_TO_COL->at(j)))
		{
			rep_c = 0;
			pr_exp = COL_TO_EXP->at(ARRPOS_TO_COL->at(j));
		}
		else
			rep_c++;

		out1 << COL_TO_EXP->at(ARRPOS_TO_COL->at(j)) << "_" << rep_c + 1 << '\t';
	}

	out1 << "---\t" << "GL_LPV\t" << "---\t";
	out2 << "---\t" << "GL_LPV\t" << "---\t";

	pr_exp = "";
	rep_c = 0;

	for(unsigned int j = 0; j < DATA->at(0).V.size(); j++)
        {
                if(pr_exp != COL_TO_EXP->at(ARRPOS_TO_COL->at(j)))
                {
                        rep_c = 0;
                        pr_exp = COL_TO_EXP->at(ARRPOS_TO_COL->at(j));
                }
                else
                        rep_c++;

                out1 << "LOC_LPV_"<< COL_TO_EXP->at(ARRPOS_TO_COL->at(j)) << "_" << rep_c + 1 << '\t';
		out2 << "LOC_LPV_"<< COL_TO_EXP->at(ARRPOS_TO_COL->at(j)) << "_" << rep_c + 1 << '\t';
        }	

	out1 << endl;
	out2 << endl;

	for(unsigned int i = 0; i < DATA->size(); i++)
	{
		out1 << DATA->at(i).gene_id << '\t' << DATA->at(i).tr_id << '\t';
		out2 << DATA->at(i).gene_id << '\t' << DATA->at(i).tr_id << '\t';

		for(unsigned int j = 0; j < DATA->at(i).v_label.size(); j++)
                {
                	out1 << DATA->at(i).v_label.at(j) << '\t';
                        out2 << DATA->at(i).v_label.at(j) << '\t';
                }

		double ttot = 0, mu = 0;

		bool gflag = false;

		for(unsigned int j = 0; j < DATA->at(i).V.size(); j++)
		{
			if(DATA->at(i).V.at(j) != 0)
				gflag = true;

			out1 << DATA->at(i).V.at(j) << '\t';
//			out2 << DATA->at(i).V.at(j) << '\t';
		}
	
		out1 << "---\t";
		out2 << "---\t";
	
		if(gflag)
		{
			double gsum = 0, gchi;
			double bi = 1.0 / (double)DATA->at(i).V.size();

			for(unsigned int j = 0; j < DATA->at(i).V.size(); j++)
				ttot += DATA->at(i).V.at(j);

			mu = ttot / (double)DATA->at(i).V.size();

			for(unsigned int j = 0; j < DATA->at(i).V.size(); j++)
				if(DATA->at(i).V.at(j) != 0)
					gsum += DATA->at(i).V.at(j) * log(DATA->at(i).V.at(j) / mu); 

			gchi = gsl_cdf_chisq_Q(2*gsum, DATA->at(i).V.size() - 1);

			if(gchi > 0)
			{
				out1 << -log10(gchi) << '\t';
				out2 << -log10(gchi) << '\t';
			}
			else
			{
				out1 << BIGLOG << '\t';
				out2 << BIGLOG << '\t';
			}
		}

		else
		{
			out1 << "0\t";
			out2 << "0\t";
		}

		out1 << "---\t";
		out2 << "---\t";

		if(gflag)
		{
			for(unsigned int j = 0; j < DATA->at(i).V.size(); j++)  //LOCAL
                	{
                        	double pv = 0, lmu = 0, ltot = DATA->at(i).V.at(j), ecount = 1, rcount = 0;
				string cexp = COL_TO_EXP->at(ARRPOS_TO_COL->at(j));

				for(unsigned int k = 0; k < DATA->at(i).V.size(); k++)     //AVOID REPLICATES
					if(COL_TO_EXP->at(ARRPOS_TO_COL->at(k)) != cexp)
					{	
						ltot += DATA->at(i).V.at(k);
						ecount++;
					}
					else
						rcount++;

				lmu = ltot/ecount;
	//			cerr << endl << "V = " << DATA->at(i).V.at(j) << '\t' << "ltot = " << ltot << '\t' << "ecount = " << ecount << '\t' << "lmu = " << lmu << endl;

                        	if(DATA->at(i).V.at(j) != 0)
                                	pv  = DATA->at(i).V.at(j) * log((DATA->at(i).V.at(j)) / lmu);

                        	if(ltot != DATA->at(i).V.at(j))
                                	pv += (ltot - DATA->at(i).V.at(j)) * log(((ltot) - DATA->at(i).V.at(j))/(((double)(DATA->at(i).V.size()) - rcount) * lmu));        

                //        	DATA->at(i).VP.push_back(gsl_cdf_chisq_Q(2*pv, 1));

				double chi = gsl_cdf_chisq_Q(2*pv, 1);

				if(DATA->at(i).V.at(j) >= lmu)
				{
					if(chi > 0)
					{
                				out1 << -log10(chi) << '\t';
						out2 << -log10(chi) << '\t';
					}
					else
					{
						out1 << BIGLOG << '\t';
						out2 << BIGLOG << '\t';
					}
				}
				else
				{
					if(chi > 0)
					{
						out1 << log10(chi) << '\t';
						out2 << log10(chi) << '\t';
					}
					else
					{
						out1 << -BIGLOG << '\t';
						out2 << -BIGLOG << '\t';
					}
				}
                	}
		}
		else
		{
			for(unsigned int j = 0; j < DATA->at(i).V.size(); j++)
			{
				out1 << "0" << '\t';
				out2 << "0" << '\t';
			}
		}

		out1 << endl;
		out2 << endl;
		
	}	

	out1.close();
	out2.close();

	return;
}

/*void calculate_pv(vector<data_row> *DATA)
{

	for(unsigned int i = 0; i < DATA->size(); i++)
	{
		double P = 0;
		double vtot = 0;

		for(unsigned int j = 0; j < DATA->at(i).V.size(); j++)
		{
			if(DATA->at(i).V.at(j) != 0)
			{
				P += DATA->at(i).V.at(j) * log(DATA->at(i).V.at(j) / DATA->at(i).mu);
				vtot += DATA->at(i).V.at(j);
			}
		}

		P *= 2;	

		DATA->at(i).pv = gsl_cdf_chisq_Q(P, DATA->at(i).V.size() - 1);

		for(unsigned int j = 0; j < DATA->at(i).V.size(); j++)
		{
			double pv = 0;

			if(DATA->at(i).V.at(j) != 0)
				pv  = DATA->at(i).V.at(j) * log((DATA->at(i).V.at(j)) / DATA->at(i).mu);

			if(vtot != DATA->at(i).V.at(j))
				pv += (vtot - DATA->at(i).V.at(j)) * log(((vtot) - DATA->at(i).V.at(j))/(((double)DATA->at(i).V.size() - 1.0)* DATA->at(i).mu));	

			DATA->at(i).VP.push_back(gsl_cdf_chisq_Q(2*pv, 1));
		}
	}


	return;
}*/

void get_data_from_file(string *data_file, vector<bool> *vb_col, vector<bool> *vb_col_lab, vector<data_row> *DATA, map<string, vector<unsigned int> > *GENE_MAP, map<unsigned int, unsigned int> *ARRPOS_TO_COL)
{
	ifstream in(data_file->c_str());

	if(!in)
        {
                cerr << "\n\nCan't find file: " << *data_file << endl;
                exit(EXIT_FAILURE);
        }

        string line;

        while(getline(in,line))
	{
		if(line.empty())
                	continue;
                if(line[0] == '#')
                	continue;

                data_row D;

		parse_row(&line, vb_col, vb_col_lab, &D, ARRPOS_TO_COL);

                if(D.valid)
		{
			if(!D.tr_id.empty() && !D.gene_id.empty())
			{
				map<string,vector<unsigned int> >::iterator mi = GENE_MAP->find(D.gene_id);

				if(mi == GENE_MAP->end())
					mi = GENE_MAP->insert(make_pair(D.gene_id, vector<unsigned int>())).first;

				mi->second.push_back(DATA->size());
			}

                        DATA->push_back(D);
		}
                else
                        cerr << endl << endl << "Invalid data found in line:\n" << line << endl << "Skipping it..." << endl; 
	}

	in.close();

	return;
}

void parse_row(string *line, vector<bool> *vb_col, vector<bool> *vb_col_lab, data_row *D, map<unsigned int, unsigned int> *ARRPOS_TO_COL)
{
	D->lrow = *line;
        D->valid = true;
//	D->mu = 0;
        vector<bool> vb_check(vb_col->size(), false);

	istringstream str(D->lrow);

	unsigned int c_count = 1;
        string buf;

	while(str >> buf)
	{
		if(c_count < vb_col->size() || c_count < vb_col_lab->size())
		{
			if(vb_col->at(c_count))
                        {
                                double v;
				vb_check.at(c_count) = true;
                                v = atof(buf.c_str());

				ARRPOS_TO_COL->insert(make_pair(D->V.size(), c_count));
				D->V.push_back(v);		

				if(buf.find_first_not_of("-01234567890eE.") != string::npos || buf.empty())
                                        D->valid = false;

	//			if(D->valid)
	//				D->mu += v;
			}

			if(vb_col_lab->at(c_count))
                                D->v_label.push_back(buf);
		}

		if(c_count == TR_COL)
			D->tr_id = buf;
		if(c_count == GENE_COL)
			D->gene_id = buf;

		c_count++;
	}

//	if(D->valid)
//		D->mu /= (double)D->V.size();

	for(unsigned int i = 0; i < vb_col->size(); i++)
        {
                if(vb_col->at(i) != vb_check.at(i))
                {
                        D->valid = false;
                        return;
                }
        }

	return;
}

void command_line_parser(int argc, char **argv)
{
        if(argc != 3)
        {
                cerr << "\nSYNTAX:\n\nRNentropy -f data_file\n\n";
                exit(EXIT_FAILURE);
        }

        for(int i = 1; i < argc; i++)
        {
                string buf = argv[i];

                if(buf == "-f")
                {
                        if(i < argc - 1)
                                DATA_FILE = argv[++i];
                        continue;
                }

     /*         else if(buf == "-c")
                {
                        if(i < argc - 1)
                                COL = argv[++i];
                        continue;
                }
                else if(buf == "-l")
                {
                        if(i < argc - 1)
                                LCOL = argv[++i];
                        continue;
                }
		else if(buf == "-tc")
		{
			if(i < argc - 1)	
				TR_COL = atoi(argv[++i]);
			continue;
		}
		else if(buf == "-gc")
		{
			 if(i < argc - 1)
				GENE_COL = atoi(argv[++i]);
			continue;
		}
		*/
                else
                {
                        cerr << "\nUnknown option: " << buf << endl << endl;
                        exit(1);
                }
        }

	return;
}

/*void col_num_parser(vector<bool> *VB_COL)
{
	for(unsigned int i = 0; i < COL.size(); i++)
	{
		if(isdigit(COL[i]))
			continue;

		else if(COL[i] == ',')
			COL[i] = ' ';

		else
		{
			cerr << "\nERROR: the correct format of the columns parameters (-c) is col_num_1,col_num_2,...,col_num_N\nExample: 1,3,6\tColumn numeration is 1 based.\n\n";
			exit(EXIT_FAILURE); 
		}
	}

	set<unsigned int> s_col;
	istringstream str(COL);
	unsigned int buf, col_num;

	while(str >> buf)
	{
		if(!buf)
		{
			cerr << "\nERROR: Column 0 not allowed.\n\n";
                        exit(EXIT_FAILURE);
		}
		
		s_col.insert(buf);
	}
	
	if(s_col.empty())
	{
		cerr << "\nERROR: Columns parameter (-c) is missing...\n\n";
		exit(EXIT_FAILURE);
	}

	col_num = *--s_col.end();	
	
	VB_COL->resize(col_num + 1, 0);

	for(unsigned int i = 1; i <= col_num; i++)
		if(s_col.find(i) != s_col.end())
			VB_COL->at(i) = true;

	cerr << endl;
        for(unsigned int i = 0; i < VB_COL->size(); i++)
                cerr << VB_COL->at(i) << ' ';

        cerr << endl;

	return;
}*/

void col_num_parser(vector<bool> *VB_COL)
{
	unsigned int col_num;

        if(S_COL.empty())
        {
                cerr << "\nERROR: No data columns...\n\n";
                exit(EXIT_FAILURE);
        }

        col_num = *--S_COL.end();

//	cerr << endl << "col_num = " << col_num << endl;

        VB_COL->resize(col_num + 1, 0);

        for(unsigned int i = 1; i <= col_num; i++)
                if(S_COL.find(i) != S_COL.end())
                        VB_COL->at(i) = true;

/*	cerr << endl;
	for(unsigned int i = 0; i < VB_COL->size(); i++)
		cerr << VB_COL->at(i) << ' ';

	cerr << endl;
*/
        return;
}

void col_lab_parser(vector<bool> *VB_COL)
{
        if(LCOL.empty())
                return;

        for(unsigned int i = 0; i < LCOL.size(); i++)
        {
                if(isdigit(LCOL[i]))
                        continue;

                else if(LCOL[i] == ',')
                        LCOL[i] = ' ';

                else
                {
                        cerr << "\nERROR: the correct format of the label columns parameters (-l) is col_num_1,col_num_2,...,col_num_N\nExample: 1,3,6\tColumn numeration is 1 based.\n\n";
                        exit(EXIT_FAILURE);
                }
        }

        set<unsigned int> s_col;
        istringstream str(LCOL);
        unsigned int buf, col_num;

        while(str >> buf)
        {
                if(!buf)
                {
                        cerr << "\nERROR: Column 0 not allowed.\n\n";
                        exit(EXIT_FAILURE);
                }

                s_col.insert(buf);
        }

        if(s_col.empty())
                return;

        col_num = *--s_col.end();

        VB_COL->resize(col_num + 1, 0);

        for(unsigned int i = 1; i <= col_num; i++)
                if(s_col.find(i) != s_col.end())
                        VB_COL->at(i) = true;

        return;
}
