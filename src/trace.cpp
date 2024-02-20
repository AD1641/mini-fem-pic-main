/*==============================================================================*
 *
 * TRACE
 *------------------------------------------------------------------------------*
 * Maintainer: Ed Higgins <ed.higgins@york.ac.uk>
 *------------------------------------------------------------------------------*
 * Version: 0.1.1, 2022-10-05
 *------------------------------------------------------------------------------*
 * This code is distributed under the MIT license.
 *==============================================================================*/

#include "trace.h"
#include <string>

#define DEBUG false

double _walltime() {
    return omp_get_wtime();
}

double Trace::calculation_time() {
    return _walltime() - trace_start;
}

Trace::Trace(std::string name_): name(std::move(name_)) {
    trace_start = _walltime();
};

void Trace::write_callstack() {
    std::cout << "Callstack:" << std::endl;
    for(const auto &call: callstack) {
        std::cout << "    " << call->name << "\t" << call->total_time << std::endl;
    }
};

void Trace::write_profile(std::string filename) {
    auto calculation_time = _walltime() - trace_start;
    std::ofstream trace_file;
    trace_file.open(filename);
    trace_file << "function_name"<< ","
               << "total_time" << ","
               << "exclusive_time" << ","
               << "exclusive_proportion" << ","
               << "number_of_calls"
               << std::endl;

    for(const auto& [func_name, func]: function_list) {
        if (!func.start_time) {
            trace_file << func_name << ","
                       << func.total_time << ","
                       << func.total_time - func.child_time << ","
                       << (func.total_time - func.child_time) / calculation_time << ","
                       << func.num_calls
                       << std::endl;
        }
    }
    trace_file << "traceEnter" << ","
               << entry_time << ","
               << entry_time << ","
               << (entry_time) / calculation_time << ","
               << 1
               << std::endl;

    trace_file << "traceExit" << ","
               << exit_time << ","
               << exit_time << ","
               << (entry_time) / calculation_time << ","
               << 1
               << std::endl;

    trace_file << "traceTotalTime" << ","
               << calculation_time << ","
               << 0 << ","
               << 0 << ","
               << 1
               << std::endl;
    trace_file.close();
};

void Trace::enter(std::string func_name) {
    if (DEBUG) std::cout << "Entering " << func_name << std::endl;

    auto start_time = _walltime();


    // if (!function_list.contains(func_name)) {
    if (function_list.find(func_name) == function_list.end()) {
        function_list.insert(std::make_pair(func_name, TraceElement()));
        function_list[func_name].name = func_name;
    }

    callstack.push_back(&function_list[func_name]);
    function_list[func_name].num_calls += 1;
    function_list[func_name].start_time = start_time;

    entry_time += (_walltime() - start_time);
};

void Trace::exit(std::string func_name) {
    if (DEBUG) std::cout << "Exiting " << func_name << std::endl;
    auto end_time = _walltime();

    auto total_time = end_time - function_list[func_name].start_time;
    function_list[func_name].total_time += total_time;
    function_list[func_name].start_time = 0;
    callstack.pop_back();
    if (!callstack.empty()) callstack.back()->child_time += total_time;
    exit_time += (_walltime() - end_time);
}

TraceCaller::TraceCaller(std::string name_) : name(std::move(name_)) {
    if (trace::enabled) {
        trace::current.enter(name);
    }
}

TraceCaller::~TraceCaller() {
    if (trace::enabled) {
        trace::current.exit(name);
    }
}

