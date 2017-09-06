#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <algorithm>
#include <vector>


using namespace std;


double distance(double *v1, double *v2);

int main (int argc, char* argv[])
{

//input file with pdb structure
std::string name(argv[1]);
string helix[10000]; 
string sheet[10000]; 

// variables declaration
string line;
int b=0; int e=0; int a=1; int ca=0; // b = begin; e = end; a = n atoms and ca = n calpha atoms
int iall=0; int icalpha=0;
double d, da; // distance
//

// open pdb file and store as pdb
ifstream pdb;
pdb.open(argv[1]);


// Get begin index and set position
int h=0; int s=0;	
while (line.substr(0,4)!="ATOM")
	{
	getline(pdb,line);
	b++;
	if(line.substr(10,3)=="DNA" || line.substr(10,3)=="RNA" || line.substr(10,12)=="CARBOHYDRATE")
	{
	cout << "\n\tNot a protein structure!\n" << endl;
	return 0;
	}
	if (line.substr(0,5)=="HELIX")
		{
		for (int m = stoi(line.substr(21,4));m<=std::stoi(line.substr(33,4));m++)
			{
			helix[h] = to_string(m)+line.substr(19,1);
			//cout << helix[h] << "-";
			h++;
			}
		//cout << "\n";
		}
	if (line.substr(0,5)=="SHEET")
		{
		for (int m=std::stoi(line.substr(22,4));m<=std::stoi(line.substr(33,4));m++)
			{
			sheet[s] = to_string(m)+line.substr(19,1);
			//cout << sheet[s] << "-";
			s++;
			}
		//cout << "\n";
		}
	}

// Get end index
e=b;
while (line.substr(0,3) != "END")
	{
	getline(pdb,line);
	e++;
	if (line.substr(0,4)=="ATOM" && line.substr(77,1)!="H")
		{
		a++;
		}
	if (line.substr(0,4)=="ATOM" && line.substr(13,2)=="CA")
		{
		ca++;
		}

	}
	
cout << "\nReading protein... PDB code: "+name+"\n" << endl;
//cout << "Begin of coordinates: " << b << endl;
//cout << "End of file: " << e << endl;
cout << "Number of non-redundant atoms in protein: " << a << endl;
cout << "Number of Calpha atoms: " << ca << endl;
cout << "\nCalculating physical interactions between all atoms..." << endl;
//////////////////////////////////

int number[a];
string at[50000];
string rt[50000];
string chain[50000];
int rn[a];
double x[a],y[a],z[a];

////// Return to the protein begin
pdb.seekg(ios::beg); // restart getline index
int k=0;
while (k<b)
{
getline(pdb,line);
k++;
}

//cout << line << endl;
int j=0;
while (line.substr(0,3) != "END")
	{
	if (line.substr(0,4)=="ATOM" && (line.substr(16,1)==" " || line.substr(16,1)=="A" ) && line.substr(77,1)!="H")
		{
		number[j] = std::stoi(line.substr(6,6));
		rn[j] = std::stoi(line.substr(22,4));
		at[j] = line.substr(13,2);
		rt[j] = line.substr(17,3);
		chain[j] = line.substr(21,1);
		x[j] = stod(line.substr(30,8));
		y[j] = stod(line.substr(39,8));
		z[j] = stod(line.substr(47,8));
		j++;
		}
	getline(pdb,line);
	}

///////
ofstream allatom;
allatom.open(name.substr(0,4)+"_allatom.txt", ios::out);
allatom << "resid\tch\tatom\tatype\trtype\t\tresid\tch\tatom\tatype\trtype\tdist (Ang)\tpair\t\tindex" << endl;
///////
string interact[100000]; //Initializing array containing interactions pairs

int t=0; int counter;
for (int w1 = 0; w1 < j; w1++)
	{
	double v1[3]={x[w1],y[w1],z[w1]};
	for (int w2 = w1+1; w2 < j; w2++)
		{
		double v2[3]={x[w2],y[w2],z[w2]};
		d = distance(v1,v2);
		if (d<=4.0 && rn[w2]-rn[w1]>3)
			{
			allatom << rn[w1] << "\t" << chain[w1] << "\t" << number[w1] << "\t"+at[w1]+"\t"+rt[w1]+"\t\t" << rn[w2] << "\t" << chain[w2] << "\t" << number[w2] << "\t"+at[w2]+"\t"+rt[w2]+"\t" << d << "\t";
			iall++;
			int p = 0;
			p = count(interact,interact+t,to_string(rn[w1])+"_"+to_string(rn[w2])+"_"+chain[w1]+"_"+chain[w2]);
			if (p == 0)
				{				
				interact[t]= to_string(rn[w1])+"_"+to_string(rn[w2])+"_"+chain[w1]+"_"+chain[w2];
				allatom << "\t" << to_string(rn[w1])+chain[w1]+"-"+to_string(rn[w2])+chain[w2] << "\t\t" << t << endl;
				t++;
				}
			if (p != 0)
				{
				allatom << "\t" << endl;
				}
			}
		}
	}

cout << "\nNumber of interactions found between all atoms in the structure: " << iall << endl;
cout << "File with distances saved as: "+name.substr(0,4)+"_allatom.txt\n" << endl;
cout << "Calculating Calpha distances for identified interactions..." << endl;


string code = interact[1000];
//cout << code << endl;
int pos = code.find("_");
string resid1 = code.substr(0,pos);
code = code.substr(pos+1);
pos = code.find("_");
string resid2 = code.substr(0,pos);
code = code.substr(pos+1);
pos = code.find("_");
string c1 = code.substr(0,pos);
code = code.substr(pos+1);
pos = code.find("_");
string c2 = code.substr(0,pos);

///////
ofstream alphatom;
alphatom.open(name.substr(0,4)+"_calpha.txt", ios::out);
alphatom << "resid\tch\trtype\tSS\t\tresid\tch\trtype\tSS\tCA-CA dist (Ang.)" << endl;
///////


double va[3]={0,0,0};
double vb[3]={0,0,0};

for (int ndx=0;ndx<t;ndx++)
	{
	string code = interact[ndx];
	//cout << code << endl;
	int pos = code.find("_");
	string resid1 = code.substr(0,pos);
	code = code.substr(pos+1);
	pos = code.find("_");
	string resid2 = code.substr(0,pos);
	code = code.substr(pos+1);
	pos = code.find("_");
	string c1 = code.substr(0,pos);
	code = code.substr(pos+1);
	pos = code.find("_");
	string c2 = code.substr(0,pos);

	for (int k=0;k<=j;k++)
		{
		if (to_string(rn[k])==resid1 && chain[k]==c1 && at[k] == "CA")
			{
			va[0]=x[k];
			va[1]=y[k];
			va[2]=z[k];
			alphatom << rn[k] << "\t" << chain[k] << "\t" << rt[k] << "\t";
			if (count(helix,helix+h,to_string(rn[k])+chain[k])>0)
				{
				alphatom << "H" << "\t";
				} 
			else if (count(sheet,sheet+s,to_string(rn[k])+chain[k])>0)
				{
				alphatom << "S" << "\t";
				} 				
			else {alphatom << "C" << "\t";}
			icalpha++;
			}
		if (to_string(rn[k])==resid2 && chain[k]==c2 && at[k] == "CA")
			{
			vb[0]=x[k];
			vb[1]=y[k];
			vb[2]=z[k];
			alphatom << "\t" << rn[k] << "\t" << chain[k] << "\t" << rt[k] << "\t";
			if (count(helix,helix+h,to_string(rn[k])+chain[k])>0)
				{
				alphatom << "H" << "\t";
				}
			else if (count(sheet,sheet+s,to_string(rn[k])+chain[k])>0)
				{
				alphatom << "S" << "\t";
				} 				
			else {alphatom << "C" << "\t";}
			}
		}
	alphatom << distance(va,vb) << endl;
	}
cout << "Number of calpha distances calculated for the structure: " << icalpha << endl;
cout << "\nFile with calpha distances saved as: "+name.substr(0,4)+"_calpha.txt\n" << endl;

allatom.close();
alphatom.close();

return 0;
}

// Calculate distance between a pair of atoms
double distance(double *v1, double *v2)
{
double dx = pow(v1[0]-v2[0],2);
double dy = pow(v1[1]-v2[1],2);
double dz = pow(v1[2]-v2[2],2);
double dist = pow(dx+dy+dz,0.5);
return dist;
}
