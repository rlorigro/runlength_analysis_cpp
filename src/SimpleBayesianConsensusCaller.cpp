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
#include "SimpleBayesianConsensusCaller.hpp"
#include "Base.hpp"

using std::runtime_error;
using std::to_string;
using std::cout;
using std::endl;
using std::max;
using Separator = boost::char_separator<char>;
using Tokenizer = boost::tokenizer<Separator>;


// Helper function
void SimpleBayesianConsensusCaller::splitAsDouble(string s, string& separators, vector<double>& tokens){
    Separator separator(separators.c_str());
    Tokenizer tok{s, separator};

    for (string token: tok) {
        tokens.push_back(stod(token));
    }
}


// Helper function
void SimpleBayesianConsensusCaller::splitAsString(string s, string& separators, vector<string>& tokens){
    Separator separator(separators.c_str());
    Tokenizer tok{s, separator};

    for (string token: tok) {
        tokens.push_back(token);
    }
}


void SimpleBayesianConsensusCaller::validateMatrixDimensions(){
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

SimpleBayesianConsensusCaller::SimpleBayesianConsensusCaller(
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



void SimpleBayesianConsensusCaller::printProbabilityMatrices(char separator) const{
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


void SimpleBayesianConsensusCaller::printPriors(char separator) const{
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


void SimpleBayesianConsensusCaller::parseName(ifstream& matrixFile, string& line){
    // Expect only one line to follow
    getline(matrixFile, line);
    configurationName = line;
}


void SimpleBayesianConsensusCaller::parsePrior(ifstream& matrixFile, string& line, vector<string>& tokens){
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


void SimpleBayesianConsensusCaller::parseLikelihood(ifstream& matrixFile, string& line, vector<string>& tokens){
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


void SimpleBayesianConsensusCaller::loadConfiguration(ifstream& matrixFile){
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


void SimpleBayesianConsensusCaller::printLogLikelihoodVector(vector<double>& logLikelihoods, size_t cutoff) const{
    size_t i = 0;
    for (auto& item: logLikelihoods){
        cout << i << " " << pow(10, item) << '\n';
        i++;
        if (i == cutoff){
            break;
        }
    }
}


void SimpleBayesianConsensusCaller::normalizeLikelihoods(vector<double>& x, double xMax) const{
    for (uint32_t i=0; i<x.size(); i++){
        x[i] = x[i]-xMax;
    }
}


void SimpleBayesianConsensusCaller::factorRepeats(
    array<std::map<uint16_t,uint16_t>,2>& factoredRepeats,
    const vector <vector <float> >& pileup_column) const{

    // Store counts for each unique observation
    for (auto& observation: pileup_column ){
        // If NOT a gap, always increment
        if (not is_gap(observation[BASE])) {
            factoredRepeats[uint16_t(observation[STRAND])][uint16_t(observation[LENGTH])]++;
        // If IS a gap only increment if "countGapsAsZeros" is true
        }else if (countGapsAsZeros){
            factoredRepeats[uint16_t(observation[STRAND])][0]++;
        }
    }
}


void SimpleBayesianConsensusCaller::factorRepeats(
    array<std::map<uint16_t,uint16_t>,2>& factoredRepeats,
    const vector <vector <float> >& pileup_column,
    uint8_t consensus_base_index) const{

    // Store counts for each unique observation
    for (auto& observation: pileup_column){
        // Ignore non consensus repeat values
        if (observation[BASE] == consensus_base_index){
            // If NOT a gap, always increment
            if (not is_gap(observation[BASE])) {
                factoredRepeats[uint16_t(observation[STRAND])][uint16_t(observation[LENGTH])]++;
            // If IS a gap only increment if "countGapsAsZeros" is true
            }else if (countGapsAsZeros){
                factoredRepeats[uint16_t(observation[STRAND])][0]++;
            }
        }
    }
}


uint16_t SimpleBayesianConsensusCaller::predictRunlength(
        const vector <vector <float> >& pileup_column,
        uint8_t consensus_base_index,
        vector<double>& logLikelihoodY) const{
    array <std::map <uint16_t,uint16_t>, 2> factoredRepeats;    // Repeats grouped by strand and length

    size_t priorIndex = -1;   // Used to determine which prior probability vector to access (AT=0 or GC=1)
    uint16_t x;               // Element of X = {x_0, x_1, ..., x_i} observed repeats
    uint16_t c;               // Number of times x_i was observed
    uint16_t y;               // Element of Y = {y_0, y_1, ..., y_j} true repeat between 0 and j=max_runlength
    double logSum;            // Product (in logspace) of P(x_i|y_j) for each i

    double yMaxLikelihood = -INF;     // Probability of most probable true repeat length
    uint16_t yMax = 0;                 // Most probable repeat length

    // Determine which index to use for this->priors
    if (index_to_base(consensus_base_index) == "A" || index_to_base(consensus_base_index) == "T"){
        priorIndex = 0;
    }
    else if (index_to_base(consensus_base_index) == "G" || index_to_base(consensus_base_index) == "C"){
        priorIndex = 1;
    }

    // Count the number of times each unique repeat was observed, to reduce redundancy in calculating log likelihoods/
    // Depending on class boolean "ignoreNonConsensusBaseRepeats" filter out observations
    if (ignoreNonConsensusBaseRepeats) {
        factorRepeats(factoredRepeats, pileup_column, consensus_base_index);
    }
    else {
        factorRepeats(factoredRepeats, pileup_column);
    }

    // Iterate all possible Y from 0 to j to calculate p(Y_j|X) where X is all observations 0 to i,
    // assuming i and j are less than maxRunlength
    for (y = 0; y <= maxOutputRunlength; y++){
        // Initialize logSum for this Y value using empirically determined priors
        logSum = priors[priorIndex][y];

        for (uint16_t strand = 0; strand <= factoredRepeats.size() - 1; strand++){
            for (auto& item: factoredRepeats[strand]) {
                x = item.first;
                c = item.second;

                // In the case that observed runlength is too large for the matrix, cap it at maxRunlength
                if (x > maxInputRunlength){
                    x = maxInputRunlength;
                }

                // Increment log likelihood for this y_j
                logSum += double(c)*probabilityMatrices[consensus_base_index][y][x];
            }
        }

        logLikelihoodY[y] = logSum;

        // Update max Y value if log likelihood is greater than previous maximum... Should do this outside this loop?
        if (logSum > yMaxLikelihood){
            yMaxLikelihood = logSum;
            yMax = y;
        }
    }

    normalizeLikelihoods(logLikelihoodY, yMaxLikelihood);

    return max(uint16_t(1), yMax);   // Don't allow zeroes...
}


uint8_t SimpleBayesianConsensusCaller::predictConsensusBase(const vector <vector <float> >& pileup_column) const{
    vector<uint32_t> baseCounts(5,0);
    uint32_t maxBaseCount = 0;
    uint8_t maxBase = 4;   // Default to gap in case coverage is empty (is this possible?)
    uint32_t key;

    // Count bases. If it's a gap increment placeholder 4 in baseCount vector
    for(auto& observation: pileup_column) {
        if (not is_gap(observation[BASE])) {
            key = observation[BASE];
            baseCounts[key]++;
        }
        else{
            key = 4;
            baseCounts[key]++;
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


vector <float> SimpleBayesianConsensusCaller::operator()(const vector <vector <float> >& coverage) const{
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

    return vector<float> {float(consensusBase), float(consensusRepeat)};
}
