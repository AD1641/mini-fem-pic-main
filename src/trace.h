/*==============================================================================*
 * TRACE
 *------------------------------------------------------------------------------*
 * Maintainer: Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-08-24
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#ifndef TRACE_H
#define TRACE_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <omp.h>

#define TRACE_ME TraceCaller _TRACE_OBJECT(__func__);
#define TRACE(x) TraceCaller _TRACE_OBJECT(x);

double _walltime();

struct TraceElement {
    std::string name;
    int num_calls = 0;
    double total_time = 0;
    double child_time = 0;
    double start_time = 0;
};

class Trace {
    private:
        double entry_time = 0;
        double exit_time = 0;
        double trace_start;
        std::string name;
        std::vector<TraceElement*> callstack;

    public:
        Trace(std::string name_);
        std::unordered_map<std::string, TraceElement> function_list;

        double calculation_time();
        void write_callstack();
        void write_profile(std::string filename);
        void enter(std::string func_name);
        void exit(std::string func_name);
};

struct TraceCaller {
    std::string name;
    TraceCaller(std::string name_);
    ~TraceCaller();

    void add(std::string name_);
};

namespace trace {
    extern int enabled;
    extern Trace current;
}

#endif /* !TRACE_H */
