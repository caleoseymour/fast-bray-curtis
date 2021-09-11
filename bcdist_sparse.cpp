// Cale Seymour
// University of Nevada, Las Vegas
// Sometime circa late 2019 (before the apocalypse)
// Fast implementation of bray curtis distance for use with a sparse OTU matrix.
// See the accessor R script for usage.

// Function bcdist:
// Given an integer matrix of 3 columns: sample id, SV id, and SV count,
// calculate bray curtis distance between samples.
// SPECIES AND SAMPLE IDS SHOULD START AT 1. DO NOT START THEM AT 0.

#include <Rcpp.h>

// Quick macro function to decompose a lower-triangle matrix into a set of 1d
// array coordinates (arrays in C/++ are fast if you know what you're doing)
#define LT_INDEX(i,j,N) (((N-1) * i - (i-1)*i/2 + (j-i))-1)

//[[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::NumericVector bcdist(const Rcpp::IntegerMatrix SparseOTUMatrix)
{
    // Identify the largest species index
    int nspecies = 0;
    for(int i = 0; i < SparseOTUMatrix.nrow(); i++)
        if ((SparseOTUMatrix(i, 0)) > nspecies)
                nspecies = SparseOTUMatrix(i, 0);
    
    // Identify the largest sample index.    
    int nsamples = 0;
    for(int i = 0; i < SparseOTUMatrix.nrow(); i++)
        if ((SparseOTUMatrix(i, 1)) > nsamples)
                nsamples = SparseOTUMatrix(i, 1);
    
    
    Rcpp::Rcout << "Output vector of length: " << (nsamples * (nsamples - 1)) / 2 << std::endl;
    
    // Debug output.
    Rcpp::Rcout << SparseOTUMatrix.nrow() << "x" << SparseOTUMatrix.ncol() <<
                " Matrix w/ " << nspecies << " species\n"
                << "and " << nsamples << " samples." << std::endl;
    
    // Allocate matrix (don't need to do lower triangle because memory is cheap)
    Rcpp::NumericVector SampleSums(nsamples);
    std::fill(SampleSums.begin(), SampleSums.end(), 0);
    
    // Stick counts in the matrix where they belong.
    Rcpp::Rcout << "Building counts..." << std::endl;
    
    unsigned short int *counts = new unsigned short int[nsamples * nspecies];
    for(int i = 0; i < nsamples*nspecies; i++)
            counts[i] = 0;
    
    Rcpp::Rcout << "... Count matrix built." << std::endl;
    
    for(int i = 0; i < SparseOTUMatrix.nrow(); i++)
        counts[(SparseOTUMatrix(i, 1)-1) * nspecies + (SparseOTUMatrix(i, 0)-1)] = static_cast<unsigned short int>(SparseOTUMatrix(i, 2)); 
    
    Rcpp::Rcout << "... Count matrix populated." << std::endl;
    
    // Sample sums
    for(int i = 0; i < SparseOTUMatrix.nrow(); i++)
        SampleSums(SparseOTUMatrix(i,1)-1) +=  static_cast<float>(SparseOTUMatrix(i, 2));

    Rcpp::Rcout << "Calculating bray curtis function." << std::endl;
    Rcpp::Rcout << "... Samples summed to get the partial denominator." << std::endl;
    
    // Allocate an array for the numerator of the equation.
    Rcpp::Rcout << "... Allocating numerator matrix" << std::endl;
    float (*numerator) = new float[static_cast<int>((nsamples * (nsamples - 1)) / 2)];
    for(int i = 0; i < (nsamples * (nsamples - 1)) / 2; i++)
        numerator[i] = 0;
    
    Rcpp::Rcout << "... Calculating the numerator of the bray curtis function." << std::endl;
    // Calculate numerator of bray curtis function
    for(int i = 0; i < nsamples - 1; i++)
        for(int j = i + 1; j < nsamples; j++)
            for(int k = 0; k < nspecies; k++)
                numerator[LT_INDEX(i,j,nsamples)] += abs(counts[(j*nspecies) + k] - counts[(i*nspecies) + k]);

    Rcpp::Rcout << "... Calculating the quotient, pairwise between samples." << std::endl;
    Rcpp::NumericVector quotient((nsamples * (nsamples - 1)) / 2);
    std::fill(quotient.begin(), quotient.end(), -1);
    for(int i = 0; i < nsamples - 1; i++)
        for(int j = i + 1; j < nsamples; j++)
        {
            //Rcpp::Rcout << i << "," << j << ":" << LT_INDEX(i,j,nsamples) << std::endl;
            quotient[LT_INDEX(i,j,nsamples)] = numerator[LT_INDEX(i,j,nsamples)] / (SampleSums(i) + SampleSums(j));
        }
    
    delete [] counts;
    delete [] numerator;
    return quotient;
}

// Fast combn script for generating pairwise combinations of indices between 1
// and n. Accepts integers as input (since it's n^2 space, you can just replace
// the names of the integers with their index in a vector of names for other
// data types).
// [[Rcpp::export]]
Rcpp::DataFrame fastcombni(Rcpp::IntegerVector inputVector)
{
    int len = inputVector.size();
    int retLen = len * (len-1) / 2;
    Rcpp::IntegerVector outputVector1(retLen);
    Rcpp::IntegerVector outputVector2(retLen);
    int start = 0;
    for (int i = 0; i < len; ++i)
    {
        for (int j = i+1; j < len; ++j)
        {
            outputVector1(start) = inputVector(i);
            outputVector2(start) = inputVector(j);
            ++start;
        }
    }
    return(Rcpp::DataFrame::create(Rcpp::Named("a") = outputVector1,
        Rcpp::Named("b") = outputVector2));
};
