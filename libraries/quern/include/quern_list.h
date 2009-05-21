#ifndef QUERN_LIST_H
#define QUERN_LIST_H

#include <algorithm>
#include <cstdlib>

// This is a light-weight singly linked list with half-decent
// allocation, suitable only for plain-old-data, and designed
// to make it easy to use robustly from plain C (e.g. no
// exceptions, basic memory management with malloc/free and
// returned error codes).

namespace Quern {

template<typename T>
struct ListEntry
{
   T data;
   ListEntry<T>* next;

   ListEntry(void) {}
   ListEntry(const T& data_, ListEntry<T>* next_=0)
      : data(data_), next(next_) {}
};

const int puddle_size=1020;

template<typename T>
struct Puddle
{
   int end;
   Puddle* next;
   ListEntry<T> array[puddle_size];
};

template<typename T>
struct Pool
{
   Puddle<T>* puddle_head;
   ListEntry<T>* free_head;

   Pool(void) : puddle_head(0), free_head(0) {}

   ~Pool(void)
   {
      Puddle<T>* p=puddle_head;
      while(p){
         Puddle<T>* p_next=p->next;
         std::free(p);
         p=p_next;
      }
      puddle_head=0;
      free_head=0;
   }

   ListEntry<T>* allocate(void)
   {
      // if we just have something free at hand, use it
      if(free_head){
         ListEntry<T>* p=free_head;
         free_head=free_head->next;
         return p;
      }
      // otherwise, if we know we need another puddle in the pool, try for it
      if(!puddle_head || puddle_head->end==puddle_size){
         Puddle<T>* new_puddle=(Puddle<T>*)std::malloc(sizeof(Puddle<T>));
         if(!new_puddle) return 0;
         new_puddle->end=0;
         new_puddle->next=puddle_head;
         puddle_head=new_puddle;
      }
      return &puddle_head->array[puddle_head->end++];
   }

   void release(ListEntry<T>* entry)
   {
      entry->next=free_head;
      free_head=entry;
   }
};

template<typename T>
struct ListIterator
{
   ListEntry<T>** handle; // pointer to a "next" or "head" pointer in the list

   ListIterator(ListEntry<T>** handle_=0)
      : handle(handle_)
   {}

   bool still_going(void) const
   {
      assert(handle);
      return *handle!=0;
   }

   void operator++(void)
   {
      assert(handle);
      assert(*handle);
      handle=&((*handle)->next);
   }

   T& operator*(void)
   {
      assert(handle);
      assert(*handle);
      return (*handle)->data;
   }

   T* operator->(void)
   {
      assert(handle);
      assert(*handle);
      return &(*handle)->data;
   }

   const T* operator->(void) const
   {
      assert(handle);
      assert(*handle);
      return &(*handle)->data;
   }

   const T& operator*(void) const
   {
      assert(handle);
      assert(*handle);
      return (*handle)->data;
   }
};

template<typename T>
struct List
{
   Pool<T>* pool;
   ListEntry<T>* head;
   int length;

   List(void) {}

   void init(Pool<T>* pool_)
   {
      pool=pool_;
      head=0;
      length=0;
   }

   ListIterator<T> begin(void)
   {
      return ListIterator<T>(&head);
   }

   bool empty(void)
   {
      return head==0;
   }

   // p ends up referring to the next element
   void erase(ListIterator<T>& p)
   {
      assert(p.handle);
      assert(*p.handle);
      ListEntry<T>* q=*p.handle;
      *p.handle=q->next;
      pool->release(q);
      --length;
   }

   T& front(void)
   {
      assert(head);
      return head->data;
   }

   const T& front(void) const
   {
      assert(head);
      return head->data;
   }

   // put the new element before the one p initially refers to;
   // if successful, p ends up referring to the new element
   bool insert(ListIterator<T>& p, const T& data)
   {
      ListEntry<T>* q=pool->allocate();
      if(!q) return false;
      q->next=*p.handle;
      q->data=data;
      *p.handle=q;
      ++length;
      return true;
   }

   void pop_front(void)
   {
      assert(head);
      ListEntry<T>* p=head;
      head=p->next;
      pool->release(p);
      --length;
   }

   bool push_front(const T& entry)
   {
      ListEntry<T>* p=pool->allocate();
      if(!p) return false;
      p->data=entry;
      p->next=head;
      head=p;
      ++length;
      return true;
   }

   int size(void)
   {
      return length;
   }

   void swap(List<T>& other)
   {
      assert(pool==other.pool);
      std::swap(head, other.head);
      std::swap(length, other.length);
   }
};

} // end namespace

#endif
