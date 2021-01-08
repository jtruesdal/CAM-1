#include <sys/time.h>
#include <sys/resource.h>

void get_maxrss(long *rss)
{
        struct rusage usage;

        getrusage(RUSAGE_SELF, &usage);

	*rss = usage.ru_maxrss;
}
