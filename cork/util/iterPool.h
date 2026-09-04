// +-------------------------------------------------------------------------
// | iterPool.h
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
#ifndef CORK_ITERPOOL_H_HEADER_HAS_BEEN_INCLUDED
#define CORK_ITERPOOL_H_HEADER_HAS_BEEN_INCLUDED

#include "prelude.h"
#include "memPool.h"
#include "parallel.h"
#include <utility>
#include <vector>
#include <functional>

template<class T>
class IterPool
{
private:
    uint numAlloced;
public:
    IterPool(int minInitBlocks=10) :
        numAlloced(0),
        block_list(nullptr),
        pool(minInitBlocks)
    {}
    IterPool(IterPool &&src)
        : numAlloced(src.numAlloced),
          block_list(src.block_list),
          pool(std::move(src.pool))
    {
        src.block_list = nullptr;
    }
    ~IterPool() {
        // run through and destruct all remaining elements
        for_each([](T* obj) {
            obj->~T();
        });
    }
    
    void clear() {
        for_each([](T* obj) {
            obj->~T();
        });
        numAlloced = 0;
        block_list = nullptr;
        pool.clear();
    }
    
    void operator=(IterPool &&src)
    {
        for_each([](T* obj) {
            obj->~T();
        });
        block_list = src.block_list;
        src.block_list = nullptr;
        pool = std::move(src.pool);
    }
    
    
private:
    struct Block {
        T       datum;
        Block   *next;
        Block   *prev;
    };
    Block *block_list;
    
public: // bulk allocation support
    // n contiguous objects, default constructed (in parallel), linked into
    // the iteration list ahead of everything already allocated.
    struct Bulk {
        Block *base = nullptr;
        uint   n    = 0;
        inline T* operator[](size_t i) const { return (T*)(base + i); }
    };
    Bulk alloc_bulk(uint n) {
        Bulk b;
        b.n = n;
        if(n == 0) return b;
        b.base = reinterpret_cast<Block*>(pool.allocContiguous((int)n));
        Block *base     = b.base;
        Block *old_head = block_list;
        cork_par::for_range((size_t)n, 8192, [&](size_t s, size_t e) {
            for(size_t i = s; i < e; ++i) {
                Block *blk = base + i;
                new ((T*)blk) T();
                blk->prev = (i == 0) ? nullptr : (base + i - 1);
                blk->next = (i + 1 < (size_t)n) ? (base + i + 1) : old_head;
            }
        });
        if(old_head) old_head->prev = base + n - 1;
        block_list  = base;
        numAlloced += n;
        return b;
    }

    // destroy everything and drop all memory
    void release() {
        for_each([](T* obj) { obj->~T(); });
        numAlloced = 0;
        block_list = nullptr;
        pool.release();
    }

    // Free the backing store without walking objects.  Only safe when every
    // live T is trivially destructible, or the caller already destroyed them.
    void release_memory() {
        numAlloced = 0;
        block_list = nullptr;
        pool.release();
    }

    // Parallel destructor over a contiguous alloc_bulk() range, then drop
    // memory.  The range must still be the entire live set.
    void release_bulk(Bulk b) {
        if (b.n) {
            cork_par::for_each_idx((size_t)b.n, 4096, [&](size_t i) {
                b[i]->~T();
            });
        }
        numAlloced = 0;
        block_list = nullptr;
        pool.release();
    }

public: // allocation/deallocation support
    T* alloc() {
        Block *new_block = pool.alloc();
        if(block_list) block_list->prev = new_block;
        new_block->next = block_list;
        new_block->prev = NULL;
        block_list = new_block;
        
        T* obj = (T*)new_block;
        new (obj) T(); // invoke default constructor when allocating
        
        numAlloced++;
        
        return obj;
    }
    void free(T* item) {
        if(item == NULL)   return;
        item->~T(); // invoke destructor before releasing
        
        numAlloced--;
        
        Block *ptr = (Block*)(item);
        if(ptr->next)   ptr->next->prev = ptr->prev;
        if(ptr->prev)   ptr->prev->next = ptr->next;
        if(ptr == block_list)   block_list = ptr->next;
        pool.free(ptr);
    }
    
public:
    template<class F>
    inline void for_each(F func) const {
        for(Block *block = block_list;
          block != NULL;
          block = block->next) {
            func((T*)(block));
        }
    }
    // snapshot of all live objects, in iteration order
    void collect(std::vector<T*> &out) const {
        out.clear();
        out.reserve(numAlloced);
        for(Block *block = block_list; block != NULL; block = block->next)
            out.push_back((T*)(block));
    }
    inline bool contains(T* tptr) const {
        for(Block *block = block_list;
          block != NULL;
          block = block->next) {
            if(tptr == (T*)(block))
                return true;
        }
        return false;
    }
    inline uint size() const {
        return numAlloced;
    }
public: // iteration support
    class iterator {
    public:
        iterator() : ptr(NULL) {}
        iterator(Block *init) : ptr(init) {}
        iterator(const iterator &cp) : ptr(cp.ptr) {}
        
        iterator& operator++() { // prefix version
            ptr = ptr->next;
            return *this;
        }
        iterator operator++(int) {
            iterator it(ptr);
            ptr = ptr->next;
            return it;
        } // postfix version
        T& operator*() {
            return ptr->datum;
        }
        T* operator->() {
            return (T*)(ptr);
        }
        bool operator==(const iterator &rhs) {
            return ptr == rhs.ptr;
        }
        bool operator!=(const iterator &rhs) {
            return ptr != rhs.ptr;
        }
    private:
        Block *ptr;
    };
    
    iterator begin() {
        return iterator(block_list);
    }
    iterator end() {
        return iterator(NULL);
    }
    
private:
    MemPool<Block> pool;
};

#endif

