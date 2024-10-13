#include "test.h"

int
main (void)
{
  void *(*nalloc) (size_t);
  void *(*nrealloc) (void *, size_t, size_t);
  void (*nfree) (void *, size_t);

  lazer_init();

  TEST_EXPECT (strcmp (LAZER_VERSION, lazer_get_version ()) == 0);
  TEST_EXPECT (LAZER_VERSION_MAJOR == lazer_get_version_major ());
  TEST_EXPECT (LAZER_VERSION_MINOR == lazer_get_version_minor ());
  TEST_EXPECT (LAZER_VERSION_PATCH == lazer_get_version_patch ());

  lazer_set_memory_functions (NULL, NULL, NULL);
  lazer_get_memory_functions (&nalloc, &nrealloc, &nfree);

  TEST_PASS ();
}
