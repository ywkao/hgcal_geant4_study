#ifndef _TOOLBOX_H_
#define _TOOLBOX_H_
#include <TString.h>

//void print_debug_info(TString varName, unsigned long int value, bool is_end);
//void print_debug_info(TString varName, double value, bool is_end);
//void print_debug_info(TString varName, int value, bool is_end);

namespace toolbox
{
    void set_string(std::ostringstream &hnamestr, TString token, int idx)
    {
        hnamestr.clear();
        hnamestr.str(token.Data());
        hnamestr<<idx;
    }

    void print_debug_info(TString varName, unsigned long int value, bool is_end=false)
    {
        if(is_end) printf("%s = %lu\n", varName.Data(), value);
        else       printf("%s = %lu, ", varName.Data(), value);
    }
    
    void print_debug_info(TString varName, uint32_t value, bool is_end=false)
    {
        if(is_end) printf("%s = %3u\n", varName.Data(), value);
        else       printf("%s = %3u, ", varName.Data(), value);
    }
    
    void print_debug_info(TString varName, double value, bool is_end=false)
    {
        if(is_end) printf("%s = %7.3f\n", varName.Data(), value);
        else       printf("%s = %7.3f, ", varName.Data(), value);
    }
    
    void print_debug_info(TString varName, int value, bool is_end=false)
    {
        if(is_end) printf("%s = %3d\n", varName.Data(), value);
        else       printf("%s = %3d, ", varName.Data(), value);
    }
}

#endif
