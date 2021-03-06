#include "Base.hpp"
#include "Pileup.hpp"
#include <cctype>

char complement_base(char base){
    base = toupper(base);

    if (base == 'A'){
        return 'T';
    }
    else if (base == 'C'){
        return 'G';
    }
    else if (base == 'G'){
        return 'C';
    }
    else if (base == 'T'){
        return 'A';
    }
    else{
        throw runtime_error("ERROR: invalid base has no complement: " + string(1, base));
    }
}


bool is_valid_base(string base){
    base[0] = toupper(base[0]);
    bool valid = false;

    if (base == "A"){
        valid = true;
    }
    else if (base == "C"){
        valid = true;
    }
    else if (base == "G"){
        valid = true;
    }
    else if (base == "T"){
        valid = true;
    }

    return valid;
}


bool is_valid_base(char base){
    base = toupper(base);
    bool valid = false;

    if (base == 'A'){
        valid = true;
    }
    else if (base == 'C'){
        valid = true;
    }
    else if (base == 'G'){
        valid = true;
    }
    else if (base == 'T'){
        valid = true;
    }

    return valid;
}


uint8_t base_to_index(string base){
    base[0] = toupper(base[0]);
    uint8_t index;

    if (base == "A"){
        index = 0;
    }
    else if (base == "C"){
        index = 1;
    }
    else if (base == "G"){
        index = 2;
    }
    else if (base == "T"){
        index = 3;
    }
    else{
        throw runtime_error("ERROR: base_to_index encountered invalid base: " + base);
    }

    return index;
}


uint8_t base_to_index(char base){
    base = toupper(base);
    uint8_t index;

    if (base == 'A'){
        index = 0;
    }
    else if (base == 'C'){
        index = 1;
    }
    else if (base == 'G'){
        index = 2;
    }
    else if (base == 'T'){
        index = 3;
    }
    else{
        throw runtime_error("ERROR: base_to_index encountered invalid base: " + string(1,base));
    }

    return index;
}


float base_to_float(char base){
    base = toupper(base);
    float index;

    if (base == 'A'){
        index = 0;
    }
    else if (base == 'C'){
        index = 1;
    }
    else if (base == 'G'){
        index = 2;
    }
    else if (base == 'T'){
        index = 3;
    }
    else if (base == 'N'){
        index = 4;
    }
    else{
        throw runtime_error("ERROR: base_to_float encountered invalid base: " + string(1,base));
    }

    return index;
}


string index_to_base(uint8_t index){
    return index_to_base_map[index];
}


string float_to_base(float index){
    return index_to_base_map[size_t(index)];
}


char float_to_base_char(float index){
    return index_to_base_char_map[size_t(index)];
}


bool is_valid_base_index(float index){
    return (index >= 0) and (index <=3);
}


bool is_valid_base_index(uint64_t index){
    return (index >= 0) and (index <=3);
}


bool is_valid_base_index(uint8_t index){
    return (index >= 0) and (index <=3);
}


bool is_gap(uint8_t index){
    return (index > 3);
}


bool is_gap(float index){
    return (index > 3);
}


bool is_empty(float index){
    return (index == Pileup::EMPTY);
}
