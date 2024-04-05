#pragma once

#include<CL/sycl.hpp>

#include "parameters.h"
#include "trace.h"
#include "maths.h"
#include "particles.h"
#include "meshes.h"
#include "FESolver.h"

#define MAX_SIZE 1000
template <typename T>
class kernel_array_class
{
private:
    T cache[MAX_SIZE] = {};
    int array_size = 0;
    
    extern SYCL_EXTERNAL bool valid ( const T value )
    {
        return value != T ( );
    }
    
public:
    
    extern SYCL_EXTERNAL bool push ( const T value )
    {
        if ( valid ( value ) && array_size < MAX_SIZE )
        {
            const int index = array_size + 1;
            cache[index] = value;
            array_size += 1;
            return true;
        }
        return false;
    }
    
    extern SYCL_EXTERNAL void clear ( )
    {
        for ( T &v : cache )
        {
            v = T ( );
        }
        array_size = 0;
    }
    
    extern SYCL_EXTERNAL T get ( int index )
    {
        return cache[index];
    }
    
    extern SYCL_EXTERNAL bool erase ( int index )
    {
        if ( valid ( cache[index] )  && array_size > 0 ) 
        {
            cache[index] = T ( );
            array_size -= 1;
            return true;
        }
        return false;
    }
    
    extern SYCL_EXTERNAL int size ( )
    {
        return array_size;
    }
};