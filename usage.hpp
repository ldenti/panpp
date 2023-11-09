#ifndef PANPP_USAGE
#define PANPP_USAGE

static const char *const VERSION = "PANPP, v0.0.1";

static const char *const MAIN_USAGE = "Usage: PANPP [index|search|version] -h";

static const char *const INDEX_USAGE_MESSAGE =
    "Usage: PANPP index FA+\n"
    "      -i       prefix for the index (default: RLCSA)\n"
    "      -v       verbose mode (default: false)\n"
    "      -@       set threads (default: 1)\n"
    "      -h       display this help and exit\n"
    "\n";

static const char *const SEARCH_USAGE_MESSAGE =
    "Usage: PANPP search INDEX QUERY_STRING\n"
    "      -v       verbose mode (default: false)\n"
    "      -h       display this help and exit\n"
    "\n";

#endif
