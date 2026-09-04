// +-------------------------------------------------------------------------
// | shortVec.h
// | 
// | Author: Gilbert Bernstein
// +-------------------------------------------------------------------------
// | COPYRIGHT:
// |    Copyright Gilbert Bernstein 2013
// |    See the included COPYRIGHT file for further details.
// |    
// |    This file is part of the Cork library.
// |
// |    Cork is free software: you can redistribute it and/or modify
// |    it under the terms of the GNU Lesser General Public License as
// |    published by the Free Software Foundation, either version 3 of
// |    the License, or (at your option) any later version.
// |
// |    Cork is distributed in the hope that it will be useful,
// |    but WITHOUT ANY WARRANTY; without even the implied warranty of
// |    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// |    GNU Lesser General Public License for more details.
// |
// |    You should have received a copy 
// |    of the GNU Lesser General Public License
// |    along with Cork.  If not, see <http://www.gnu.org/licenses/>.
// +-------------------------------------------------------------------------
#pragma once

#include <algorithm>
#include <new>
#include <cstring>

#include "prelude.h"

// Small vector with LEN entries of inline storage.
//
// This replaces the original implementation, which drew its small blocks
// from a *shared static MemPool*; that pool was not thread safe, which made
// it impossible to construct or destroy any topology object from more than
// one thread.  Inline storage is thread safe, removes a pointer chase on
// every access, and keeps the exact same API.

template<class T, uint LEN>
class ShortVec
{
public: // constructor/destructor
    ShortVec(uint size = 0);
    ShortVec(uint size, const T &fill_val);
    ShortVec(const ShortVec<T,LEN> &cp);
    ~ShortVec();
    
    ShortVec<T,LEN>& operator=(const ShortVec<T,LEN> &vec);
    
public: // index accessors
    inline       T& operator[](uint i)       { return data[i]; }
    inline const T& operator[](uint i) const { return data[i]; }
    
public: // iterators
    typedef       T* iterator;
    typedef const T* const_iterator;
    iterator        begin()       { return data; }
    const_iterator  begin() const { return data; }
    iterator        end()         { return data + user_size; }
    const_iterator  end()   const { return data + user_size; }
    
public: // inspectors
    inline uint size() const { return user_size; }
    
public: // modifiers
    void resize(uint newsize);
    void push_back(const T &datum);
    void erase(const T &val); // erase if it can be found
    
private: // helper functions
    inline T* inlineData() { return reinterpret_cast<T*>(inline_buf); }
    T*   allocData(uint space, uint &allocated);
    void deallocData(T* data_ptr, uint allocated);
    
    void constructRange(T* array, int begin, int end);
    void copyConstructRange(const T* src, T* dest, int begin, int end);
    void destructRange(T* array, int begin, int end);
    
    // resize, manage allocation/deallocation,
    // but not construction/destruction
    void resizeHelper(uint newsize);
    
private: // instance data
    uint user_size;     // actual number of entries from client perspective
    uint internal_size; // number of entries allocated (>= LEN)
    T* data;
    alignas(T) byte inline_buf[sizeof(T)*LEN];
};

template<class T, uint LEN> inline
T* ShortVec<T,LEN>::allocData(uint space, uint &allocated)
{
    if(space <= LEN) {
        allocated = LEN;
        return inlineData();
    }
    allocated = space;
    return reinterpret_cast<T*>(new byte[sizeof(T)*space]);
}
template<class T, uint LEN> inline
void ShortVec<T,LEN>::deallocData(T* data_ptr, uint allocated)
{
    if(allocated > LEN)
        delete[] reinterpret_cast<byte*>(data_ptr);
}

template<class T, uint LEN> inline
void ShortVec<T,LEN>::constructRange(T* array, int begin, int end)
{
    for(int i=begin; i<end; i++)
        new (&(array[i])) T();
}
template<class T, uint LEN> inline
void ShortVec<T,LEN>::copyConstructRange(const T* src, T* dest, int begin, int end)
{
    for(int i=begin; i<end; i++)
        new (&(dest[i])) T(src[i]);
}
template<class T, uint LEN> inline
void ShortVec<T,LEN>::destructRange(T* array, int begin, int end)
{
    for(int i=begin; i<end; i++)
        (&(array[i]))->~T();
}

template<class T, uint LEN> inline
void ShortVec<T,LEN>::resizeHelper(uint newsize) {
    if(newsize > internal_size) { // we need more space!
        uint new_space;
        T *new_data = allocData(std::max(newsize, internal_size*2), new_space);
        copyConstructRange(data, new_data, 0, user_size);
        destructRange(data, 0, user_size);
        deallocData(data, internal_size);
        data = new_data;
        internal_size = new_space;
    }
    user_size = newsize;
}

template<class T, uint LEN> inline
ShortVec<T,LEN>::ShortVec(uint size) : user_size(size)
{
    data = allocData(user_size, internal_size);
    constructRange(data, 0, user_size);
}
template<class T, uint LEN> inline
ShortVec<T,LEN>::ShortVec(uint size, const T &fill_val) : user_size(size)
{
    data = allocData(user_size, internal_size);
    for(uint i=0; i<user_size; i++)
        new (&data[i]) T(fill_val);
}
template<class T, uint LEN> inline
ShortVec<T,LEN>::ShortVec(const ShortVec<T,LEN> &cp) : user_size(cp.user_size)
{
    data = allocData(user_size, internal_size);
    copyConstructRange(cp.data, data, 0, user_size);
}
template<class T, uint LEN> inline
ShortVec<T,LEN>::~ShortVec()
{
    destructRange(data, 0, user_size);
    deallocData(data, internal_size);
}

template<class T, uint LEN> inline
ShortVec<T,LEN>& ShortVec<T,LEN>::operator=(const ShortVec<T,LEN> &vec)
{
    if(this == &vec) return *this;
    uint old_size = user_size;
    resizeHelper(vec.user_size);
    for(uint i=0; i<std::min(vec.user_size, old_size); i++)
        data[i] = vec.data[i];
    if(vec.user_size > old_size)
        copyConstructRange(vec.data, data, old_size, vec.user_size);
    if(vec.user_size < old_size)
        destructRange(data, vec.user_size, old_size);
    return *this;
}

template<class T, uint LEN> inline
void ShortVec<T,LEN>::resize(uint newsize) {
    uint oldsize = user_size;
    resizeHelper(newsize);
    if(oldsize < newsize)
        constructRange(data, oldsize, newsize);
    if(newsize < oldsize)
        destructRange(data, newsize, oldsize);
}

template<class T, uint LEN> inline
void ShortVec<T,LEN>::push_back(const T &datum)
{
    uint i = user_size;
    resizeHelper(user_size+1); // make room
    new (&(data[i])) T(datum);
}

template<class T, uint LEN> inline
void ShortVec<T,LEN>::erase(const T &val)
{
    for(uint i=0; i<user_size; i++) {
        if(data[i] == val) {
            std::swap(data[i], data[user_size-1]);
            resize(user_size-1);
            break;
        }
    }
}
