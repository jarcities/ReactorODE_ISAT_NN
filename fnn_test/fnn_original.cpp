#include <cmath>
#include <iostream>
#include <cstdlib> 

namespace Gl
{

    int nLayers = 5;
    int nNeurons = 10;
    int nx = 11;

    int ia[100];
    int ib[100];
    int n1[100];
    int n2[100];

    double A[1000000];
    double b[10000];

    void initfnn()
    {

        int i1 = 0;
        int i2 = 0;

        char file1[50];
        char file2[50];

        for (int ll = 1; ll <= nLayers; ll++)
        {

            ia[ll - 1] = i1;
            ib[ll - 1] = i2;

            sprintf(file1, "./A%d.csv", ll);
            sprintf(file2, "./B%d.csv", ll);

            n1[ll - 1] = nNeurons;
            n2[ll - 1] = nNeurons;
            if (ll == 1)
            {
                n1[ll - 1] = nx;
            }
            if (ll == nLayers)
            {
                n2[ll - 1] = nx;
            }

            FILE *pFile;
            float a;
            pFile = fopen(file1, "r+");
            for (int ii = 0; ii < n1[ll - 1] * n2[ll - 1]; ii++)
            {
                fscanf(pFile, "%f", &a);
                A[i1] = a;
                i1 = i1 + 1;
            }
            fclose(pFile);

            pFile = fopen(file2, "r+");
            for (int ii = 0; ii < n2[ll - 1]; ii++)
            {
                fscanf(pFile, "%f", &a);
                b[i2] = a;
                i2 = i2 + 1;
            }
            fclose(pFile);
        }

        // std::cerr<<A[ia[2]+3]<<std::endl;
        // std::cerr<<b[ib[1]+4]<<std::endl;
    }

}

using namespace Gl;

double fAct(double x)
{

    return x * tanh(log(1.0 + exp(x)));
}

// Nueral network function
void myfnn(int &nx, double x[], double fnn[])
{

    static int bbbb;

    double x1[100];
    double x2[100];

    if (bbbb != 7777)
    {
        Gl::initfnn();
        bbbb = 7777;
    }

    for (int ii = 0; ii < n1[0]; ii++)
    {
        x1[ii] = x[ii]; // initialize the input
    }

    // go through each layer of nodes
    for (int ll = 0; ll < nLayers; ll++)
    {

        // go through each node
        for (int kk = 0; kk < n2[ll]; kk++)
        {
            x2[kk] = 0.0;

            // matrix multiplication with weights 
            for (int jj = 0; jj < n1[ll]; jj++)
            {
                // x2[kk] += A[ ia[ll] + jj + (kk-1)*n1[ll] ]*x1[jj];
                x2[kk] += A[ia[ll] + kk + (jj - 1) * n2[ll]] * x1[jj];
            }

            // add bias
            x2[kk] += b[ib[ll] + kk];

            // at end of each layer, apply activation function
            if (ll < nLayers - 1)
            {
                x2[kk] = fAct(x2[kk]);
            }
        }

        // transfer outputs to inputs for next layer
        for (int kk = 0; kk < n2[ll]; kk++)
        {
            x1[kk] = x2[kk];
        }
    }

    // stores finished output in fnn
    for (int kk = 0; kk < nx; kk++)
    {
        fnn[kk] = x2[kk];
    }
}

int main()
{
    std::cout << "output of fnn" << std::endl;

    const int inputSize = Gl::nx; //nx = 11
    double x[inputSize]; //x is 1:1:11
    double fnn_out[inputSize];

    for (int i = 0; i < inputSize; i++) {
        x[i] = static_cast<double>(i);
    }

    myfnn(Gl::nx, x, fnn_out);

    for (int i = 0; i < inputSize; i++) {
        std::cout << "fnn_out[" << i << "] = " << fnn_out[i] << std::endl;
    }
}
