#include "lazer.h"

/* updated on release. XXX */
#define VERSION_MAJOR 0
#define VERSION_MINOR 1
#define VERSION_PATCH 0
#define VERSION "0.1.0"

unsigned int
lazer_get_version_major (void)
{
  return VERSION_MAJOR;
}

unsigned int
lazer_get_version_minor (void)
{
  return VERSION_MINOR;
}

unsigned int
lazer_get_version_patch (void)
{
  return VERSION_PATCH;
}

const char *
lazer_get_version (void)
{
  return VERSION;
}
