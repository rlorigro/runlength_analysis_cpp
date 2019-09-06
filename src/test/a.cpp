
#include <fstream>
#include <iostream>


using namespace std;


int main(){

    ifstream f("/home/ryan/code/runlength_analysis_cpp/output/test_CompressedRunnieWriter.rq");
    f.seekg (0, f.end);
    int length = f.tellg();
    cout << length << endl;

    return 0;
}
