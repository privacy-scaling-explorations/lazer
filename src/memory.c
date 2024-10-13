#include "memory.h"
#include "lazer.h"

#include <stddef.h>
#include <stdlib.h>

void
lazer_set_memory_functions (void *(*nalloc) (size_t),
                            void *(*nrealloc) (void *, size_t, size_t),
                            void (*nfree) (void *, size_t))
{
  if (nalloc == NULL)
    nalloc = alloc_default;
  if (nrealloc == NULL)
    nrealloc = realloc_default;
  if (nfree == NULL)
    nfree = free_default;

  _alloc = nalloc;
  _realloc = nrealloc;
  _free = nfree;
}

void
lazer_get_memory_functions (void *(**nalloc) (size_t),
                            void *(**nrealloc) (void *, size_t, size_t),
                            void (**nfree) (void *, size_t))
{
  if (nalloc != NULL)
    *nalloc = _alloc;
  if (nrealloc != NULL)
    *nrealloc = _realloc;
  if (nfree != NULL)
    *nfree = _free;
}

static void *
alloc_default (size_t len)
{
  void *mem;

  ASSERT_ERR (len > 0);

  mem = malloc (len);
  ERR (mem == NULL, "malloc failed (size %llu).", (unsigned long long)len);
  return mem;
}

static void *
realloc_default (void *mem, size_t olen, size_t nlen)
{
  ASSERT_ERR (mem != NULL);
  ASSERT_ERR (olen > 0);
  ASSERT_ERR (nlen > 0);

  mem = realloc (mem, nlen);
  ERR (mem == NULL, "realloc failed (old size %llu, new size %llu).",
       (unsigned long long)olen, (unsigned long long)nlen);
  return mem;
}

static void
free_default (void *mem, UNUSED size_t len)
{
  free (mem);
}
