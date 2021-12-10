#include <map>
#include <set>
#include <string>


class wvStringSet
{
public:
  std::set<std::string> string_set;
};


extern "C" void wv_stringSetOpen(wvStringSet **set)
{
    *set = new wvStringSet;
    (*set)->string_set.clear();
}

extern "C" void wv_stringSetClose(wvStringSet *set)
{
    set->string_set.clear();
    delete set;
}

extern "C" void wv_stringSetReset(wvStringSet *set)
{
    set->string_set.clear();
}

extern "C" int wv_stringSetContains(wvStringSet *set, const char *str)
{
    if (set->string_set.find(str) != set->string_set.end()) return 1;
    return 0;
}

extern "C" int wv_stringSetAdd(wvStringSet *set, const char *str)
{
    if (set->string_set.find(str) != set->string_set.end()) return 0;
    set->string_set.insert(str);
    return 1;
}

extern "C" int wv_stringSetDelete(wvStringSet *set, const char *str)
{
    return (int) set->string_set.erase(str);
}

extern "C" int wv_stringSetSize(wvStringSet *set)
{
    return (int) set->string_set.size();
}
