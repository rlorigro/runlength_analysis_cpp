#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <array>
#include <cmath>
#include <map>
#include "SimpleBayesianRunnieConsensusCaller.hpp"
#include "DiscreteWeibull.hpp"
#include "Miscellaneous.hpp"
#include "Base.hpp"

using std::runtime_error;
using std::to_string;
using std::cout;
using std::endl;
using std::max;
using std::log10;
using Separator = boost::char_separator<char>;
using Tokenizer = boost::tokenizer<Separator>;


// Helper function
void SimpleBayesianRunnieConsensusCaller::splitAsDouble(string s, string& separators, vector<double>& tokens){
    Separator separator(separators.c_str());
    Tokenizer tok{s, separator};

    for (string token: tok) {
        tokens.push_back(stod(token));
    }
}


// Helper function
void SimpleBayesianRunnieConsensusCaller::splitAsString(string s, string& separators, vector<string>& tokens){
    Separator separator(separators.c_str());
    Tokenizer tok{s, separator};

    for (string token: tok) {
        tokens.push_back(token);
    }
}


void SimpleBayesianRunnieConsensusCaller::validateMatrixDimensions(){
    size_t ySize = 0;
    size_t xSize = 0;

    for (size_t baseIndex=0; baseIndex < probabilityMatrices.size(); baseIndex++){
        if (baseIndex > 0){
            if (ySize != probabilityMatrices[baseIndex].size()){
                string base = index_to_base(uint8_t(baseIndex));
                throw runtime_error("ERROR: matrix size conflict detected. Matrix " + base +
                " does not match previous base.");
            }
        }

        ySize = probabilityMatrices[baseIndex].size();

        for (size_t yIndex=0; yIndex < ySize; yIndex++) {
            if (yIndex > 0) {
                if (xSize != probabilityMatrices[baseIndex][yIndex].size()) {
                    string base = index_to_base(uint8_t(baseIndex));
                    throw runtime_error("ERROR: matrix row size conflict in matrix " + base +
                    " at row " + to_string(yIndex));
                }
            }

            xSize = probabilityMatrices[baseIndex][yIndex].size();
        }
    }

    if (priors[0].size() != priors[1].size()){
        throw runtime_error("ERROR: prior probability vector sizes do not match.");
    }

    if (probabilityMatrices[0].size() != priors[0].size()){
        throw runtime_error("ERROR: prior probability vector size (" + to_string(priors[0].size()) +
        ") does not match y (true) dimension size (" + to_string(probabilityMatrices[0].size()) +
        ") of likelihood matrix.");
    }
}


// The constructor string can be either:
// - A name identifying one of the built-in configurations.
// - A path to a configuration file.

SimpleBayesianRunnieConsensusCaller::SimpleBayesianRunnieConsensusCaller(
    path& config_path){
    ignoreNonConsensusBaseRepeats = true;
    predictGapRunlengths = false;
    countGapsAsZeros = false;

    // If it was not a built-in name,
    // interpret the constructor string as a path to
    // a configuration file.
    ifstream config_file(config_path);
    loadConfiguration(config_file);

    validateMatrixDimensions();

    maxInputRunlength = uint16_t(probabilityMatrices[0][0].size() - 1);
    maxOutputRunlength = uint16_t(probabilityMatrices[0].size() - 1);

    cout << "Bayesian consensus caller configuration name is " <<
        configurationName << endl;
}



void SimpleBayesianRunnieConsensusCaller::printProbabilityMatrices(char separator) const{
    const uint32_t length = uint(probabilityMatrices[0].size());
    uint32_t nBases = 4;

    for (uint32_t b=0; b<nBases; b++){
        cout << '>' << index_to_base(uint8_t(b)) << " " << probabilityMatrices[b].size() << '\n';

        for (uint32_t i=0; i<length; i++){
            for (uint32_t j=0; j<length; j++){
                // Print with exactly 9 decimal values
                printf("%.9f",probabilityMatrices[b][i][j]);
                if (j != length-1){
                    cout << separator;
                }
            }
            cout << '\n';
        }
        if (b != nBases-1){
            cout << '\n';
        }
    }
}


void SimpleBayesianRunnieConsensusCaller::printPriors(char separator) const{
    const uint32_t length = uint(priors[0].size());
    uint32_t nBases = 2;

    for (uint32_t b=0; b<nBases; b++){
        cout << '>' << index_to_base(uint8_t(b)) << " " << priors[b].size() << '\n';

        for (uint32_t i=0; i<length; i++){
            printf("%d %.9f",int(i), priors[b][i]);
            if (i != length-1){
                cout << separator;
            }
        }
        if (b != nBases-1){
            cout << '\n';
        }
    }
}


void SimpleBayesianRunnieConsensusCaller::parseName(ifstream& matrixFile, string& line){
    // Expect only one line to follow
    getline(matrixFile, line);
    configurationName = line;
}


void SimpleBayesianRunnieConsensusCaller::parsePrior(ifstream& matrixFile, string& line, vector<string>& tokens){
    // Expect only one line to follow
    getline(matrixFile, line);

    // Initialize empty vector to fill with tokens from csv
    vector<double> row;

    // Assume csv format
    string separators = ",";
    splitAsDouble(line, separators, row);

    // Two prior distributions exist. One for AT and one for GC, since observed reads are bidirectional
    if (tokens[0] == "AT"){
        priors[0] = row;
    }
    else if (tokens[0] == "GC"){
        priors[1] = row;
    }
}


void SimpleBayesianRunnieConsensusCaller::parseLikelihood(ifstream& matrixFile, string& line, vector<string>& tokens){
    char base;
    uint32_t baseIndex = 0;
    string separators = ",";

    // Expect many lines (usually 51)
    while (getline(matrixFile, line)){

        // Stop iterating lines when blank line is encountered
        if (line.empty()){
            break;
        }

        // Initialize empty vector to fill with tokens from csv
        vector<double> row;

        // Assume csv format
        splitAsDouble(line, separators, row);

        base = tokens[0][0];
        baseIndex = uint32_t(base_to_index(base));

        probabilityMatrices[baseIndex].push_back(row);
    }
}


void SimpleBayesianRunnieConsensusCaller::loadConfiguration(ifstream& matrixFile){
    string line;
    string separators = " ";

    while (getline(matrixFile, line)){

        // Header line (labeled via fasta-like headers)
        if (line[0] == '>'){
            vector<string> tokens;

            // Store the header
            line = line.substr(1, line.size()-1);
            splitAsString(line, separators, tokens);

            if (tokens[0] == "Name"){
                parseName(matrixFile, line);

            }else if (tokens[1] == "prior"){
                parsePrior(matrixFile, line, tokens);

            }else if (tokens[1] == "likelihood"){
                parseLikelihood(matrixFile, line, tokens);
            }
        }
    }
}


void SimpleBayesianRunnieConsensusCaller::printLogLikelihoodVector(vector<double>& logLikelihoods, size_t cutoff) const{
    size_t i = 0;
    for (auto& item: logLikelihoods){
        cout << i << " " << pow(10, item) << '\n';
        i++;
        if (i == cutoff){
            break;
        }
    }
}


void SimpleBayesianRunnieConsensusCaller::normalizeLikelihoods(vector<double>& x, double xMax) const{
    for (uint32_t i=0; i<x.size(); i++){
        x[i] = x[i]-xMax;
    }
}


uint16_t SimpleBayesianRunnieConsensusCaller::predictRunlength(
        const vector <vector <float> >& pileup_column,
        uint8_t consensus_base_index,
        vector<double>& logLikelihoodY) const{

    vector<double> x(this->maxOutputRunlength+1);         // The distribution emitted by the basecaller for each pair of discrete weibull parameters
    size_t priorIndex = -1;         // Used to determine which prior probability vector to access (AT=0 or GC=1)
    uint8_t base_index;             // Base of the observation
    uint8_t forward_base_index;     // Base of the observation
    double scale;                   // The weibull parameter
    double shape;                   // The weibull parameter
    uint16_t y;                     // Element of Y = {y_0, y_1, ..., y_j} true repeat between 0 and j=max_runlength
    double logSum;                  // Product (in logspace) of P(x_i|y_j) for each i
    double log_likelihood_k;         // The kth log likelihoood in the weibull distribution of length k
    double integral_likelihood;    // Sum of P(x_i_k|y_j) for each k in 0-50 in the weibull distribution

    double yMaxLikelihood = -INF;   // Probability of most probable true repeat length
    uint16_t yMax = 0;              // Most probable repeat length

    // Determine which index to use for this->priors
    if (index_to_base(consensus_base_index) == "A" || index_to_base(consensus_base_index) == "T"){
        priorIndex = 0;
    }
    else if (index_to_base(consensus_base_index) == "G" || index_to_base(consensus_base_index) == "C"){
        priorIndex = 1;
    }

    vector<double> avg_distribution(this->maxOutputRunlength+1);
    string debug_string;

    // Iterate all possible Y from 0 to j to calculate p(Y_j|X) where X is all observations 0 to i,
    // assuming i and j are less than maxRunlength
    for (y = 0; y <= maxOutputRunlength; y++){
        // Initialize logSum for this Y value using empirically determined priors
        logSum = priors[priorIndex][y];

        for (auto& observation: pileup_column) {
            base_index = uint8_t(observation[BASE]);

            if (observation[REVERSAL] > 0){
                forward_base_index = 3 - base_index;
            }
            else{
                forward_base_index = base_index;
            }

            if (ignoreNonConsensusBaseRepeats and (forward_base_index != consensus_base_index)){
                continue;
            }

            scale = observation[SCALE];
            shape = observation[SHAPE];

            if (y==0) {
                cout << int(base_index) << " " << int(consensus_base_index) << " " << scale << " " << shape << '\n';
            }

            evaluate_discrete_weibull(x, scale, shape);

            integral_likelihood = -INF;

            // Increment log likelihood for this y_j
            for (size_t x_index=0; x_index < this->probabilityMatrices[consensus_base_index][y].size(); x_index++) {
                if (y==0) {
                    avg_distribution[x_index] += x[x_index];
                }

                if (x[x_index] == 0){
                    continue;
                }

                log_likelihood_k = log10(x[x_index]) + probabilityMatrices[consensus_base_index][y][x_index];
                integral_likelihood = log10_sum_exp(integral_likelihood, log_likelihood_k);

//                cout << x_index << " " << y << " " << x.size() << " " << probabilityMatrices.size() << " " << probabilityMatrices[y].size() << '\n';
//                debug_string += index_to_base(forward_base_index) + " " + to_string(logSum) + " " + to_string(y) + " " + to_string(x_index) + " " +  to_string(x[x_index]) + " " + to_string(log10(x[x_index])) + " " + to_string(probabilityMatrices[consensus_base_index][y][x_index]) + '\n';
            }

            logSum += integral_likelihood;
        }

        logLikelihoodY[y] = logSum;

        // Update max Y value if log likelihood is greater than previous maximum... Should do this outside this loop?
        if (logSum > yMaxLikelihood){
            yMaxLikelihood = logSum;
            yMax = y;
        }
    }

    normalizeLikelihoods(logLikelihoodY, yMaxLikelihood);

//    if (yMax == 1) {
    print_distribution(avg_distribution);
    printLogLikelihoodVector(logLikelihoodY, 15);
    cout << yMax << "\n\n";
//    }

//    if (yMax > 1){
//        cout << debug_string;
//    }

    return max(uint16_t(1), yMax);   // Don't allow zeroes...
}


uint8_t SimpleBayesianRunnieConsensusCaller::predictConsensusBase(const vector <vector <float> >& pileup_column) const{
    vector<uint32_t> baseCounts(5,0);
    uint32_t maxBaseCount = 0;
    uint8_t maxBase = 4;   // Default to gap in case coverage is empty (is this possible?)
    size_t base_index;

    // Count bases. If it's a gap increment placeholder 4 in baseCount vector
    for(auto& observation: pileup_column) {
        if (is_empty(observation[BASE])){
            continue;
        }
        if (not is_gap(observation[BASE])) {
            base_index = size_t(observation[BASE]);

            baseCounts[base_index]++;
        }
        else{
            base_index = 4;
            baseCounts[base_index]++;
        }
    }

    // Determine most represented base (consensus)
    for (uint32_t i=0; i<5; i++){
        if (baseCounts[i] > maxBaseCount){
            maxBaseCount = baseCounts[i];
            maxBase = uint8_t(i);
        }
    }

    return maxBase;
}


void SimpleBayesianRunnieConsensusCaller::operator()(const vector <vector <float> >& coverage, vector <float>& consensus) const{
    consensus = {};
    uint8_t consensusBase;
    uint16_t consensusRepeat;
    vector<double> logLikelihoods(u_long(maxOutputRunlength+1), -INF);    // initialize as zeros in log space

    consensusBase = predictConsensusBase(coverage);

    if (predictGapRunlengths) {
        // Predict all run lengths regardless of whether consensus base is a gap
        consensusRepeat = predictRunlength(coverage, consensusBase, logLikelihoods);
    }
    else{
        if (not is_gap(consensusBase)) {
            // Consensus is NOT a gap character, and the configuration forbids predicting gaps
            consensusRepeat = predictRunlength(coverage, consensusBase, logLikelihoods);
        }
        else{
            // Consensus IS a gap character, and the configuration forbids predicting gaps
            consensusRepeat = 0;
        }
    }

    consensus.emplace_back(float(consensusBase));
    consensus.emplace_back(float(consensusRepeat));
}
