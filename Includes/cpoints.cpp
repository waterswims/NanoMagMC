#include "functions.h"
#include <iostream>
#include <cmath>
#include <sstream>
#include <vector>
#include <cstdlib>


void print_cval(fstream &stream, string fname, double cval)
{
    stream.open(fname.c_str(), fstream::out | fstream::app);
    stream << cval << endl;
    stream.close();
}

void print_cval(fstream &stream, string fname, int cval)
{
    stream.open(fname.c_str(), fstream::out | fstream::app);
    stream << cval << endl;
    stream.close();
}

void print_clist(fstream &stream, string fname, double clist[], int l)
{
    stream.open(fname.c_str(), fstream::out | fstream::app);
    for(int i = 0; i < l; i++)
    {
        stream << clist[i] << " ";
    }
    stream << endl;
    stream.close();
}

int read_cval(fstream &stream, string fname, int cvals[])
{
    stream.open(fname.c_str(), fstream::in);
    int curr, i;
    for(i = 0; stream >> curr; i++)
    {
        cvals[i] = curr;
    }
    stream.close();
    return i;
}

int read_clist(fstream &stream, string fname, double clist[][100])
{
    stream.open(fname.c_str(), fstream::in);
    int i;
    string line;
    vector<string> vals;

    for(i=0; getline(stream, line); i++)
    {
        char* in = new char[line.size()+1];
        memcpy(in, line.c_str(), line.size()+1);
        char* pch;
        pch = strtok(in, " ");
        for (int j = 0; pch != NULL; j++)
        {
            clist[i][j] = atof(pch);
            pch = strtok(NULL, " ");
        }
    }
    stream.close();
    return i;
}

string cpointname(string prefix, int rank, double size, int dist, char latt, char ham, double field)
{
    stringstream stream;
    stream << "Checkpoints/" << prefix << "_" << dist << size << latt << ham << field << "_" << rank << ".cp";
    string out;
    stream >> out;
    return out;
}

bool exists(string &name)
{
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}
