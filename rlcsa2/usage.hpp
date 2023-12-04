#ifndef PANPP_USAGE
#define PANPP_USAGE

static const char *const VERSION = "PANPP, v0.0.1";

static const char *const MAIN_USAGE = "Usage: PANPP [index|search|version] -h";

static const char *const INDEX_USAGE_MESSAGE =
    "Usage: PANPP index FA+\n"
    "      -i <STR>   prefix for the index (default: RLCSA)\n"
    "      -@ <INT>   set threads (default: 1)\n"
    "      -v         verbose mode\n"
    "      -h         display this help and exit\n"
    "\n";

static const char *const SEARCH_USAGE_MESSAGE =
    "Usage: PANPP search INDEX FX\n"
    "      -f <INT>   size of flanking region on both side of a specific "
    "string (default: 0)\n"
    "      -b <INT>   batch size (default: 10000)\n"
    "      -a         do not assemble overlapping specific strings\n"
    "      -x         output in FASTX instead of TSV\n"
    "      -@ <INT>   set threads (default: 1)\n"
    "      -v         verbose mode\n"
    "      -h         display this help and exit\n"
    "\n";

#endif
